/*
  Copyright (C) 2024 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include </home/bbpdneu1/software/aspect/aspect/vp_melt_plugin/co2_melt.h>
#include <aspect/utilities.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {


      template <int dim>
      Co2Melt<dim>::Co2Melt()
        = default;


      template <int dim>
      double
      Co2Melt<dim>::
      reference_darcy_coefficient () const
      {
        // 0.01 = 1% melt
        return reference_permeability * Utilities::fixed_power<3>(0.01) / viscosity_fluid;
      }

      template <int dim>
      void
      Co2Melt<dim>::
      calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                      typename Interface<dim>::MaterialModelOutputs &out) const
      {
          ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

          for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
            {
              // Calculate new melt and compositional values. At the moment this assumes there is always
              // 3 components.

              // Define component dependent parameters 
              std::vector<double> C_bar (n_components);
              std::vector<double> composition(this->n_compositional_fields());
              const double temperature = in.temperature[q];

              // Set pressure to surface pressure if it is below zero.
              /*const double pressure    = this->get_adiabatic_conditions().pressure(in.position[q]) > 0
                                         ?
                                         this->get_adiabatic_conditions().pressure(in.position[q])
                                         :
                                         101325;*/

              double pressure = in.pressure[q] > 0
                            ? 
                            in.pressure[q]
                            :
                            101325;
              
              pressure = std::max(pressure, this->get_adiabatic_conditions().pressure(in.position[q])*0.5);

              // Get the compositional values, and limit between 0 and 1.
              double morb_cl =  std::max(0.0, std::min(in.composition[q][mcl_idx],1.0));
              double morb_cs =  std::max(0.0, std::min(in.composition[q][mcs_idx],1.0));
              double cmorb_cl =  std::max(0.0, std::min(in.composition[q][ccl_idx],1.0));
              double cmorb_cs =  std::max(0.0, std::min(in.composition[q][ccs_idx],1.0));
              double Fvol_old =  std::max(0.0, std::min(in.composition[q][melt_idx],1.0));                                                              

              // Calculate dunite and order liquid and solid components
              double dunite = 1 - morb_cs - cmorb_cs;
              double dunite_l = 1 - morb_cl - cmorb_cl;
              std::vector<double> c_s = {dunite, morb_cs, cmorb_cs};
              std::vector<double> c_l = {dunite_l, morb_cl, cmorb_cl};

              // We track the volume of melt, convert to mass here.
              double avg_rho = Fvol_old*rho_l + (1 - Fvol_old)*rho_s;
              double Fmass_old = Fvol_old*avg_rho/rho_l;

              // Now that things are ordered, find the bulk composition for each component.
              for (unsigned int i=0; i<n_components; ++i)
                C_bar[i] = Fmass_old*c_l[i] + (1-Fmass_old)*c_s[i];

              const double T_solidus = T_solidus_liquidus(pressure, C_bar, true);
              const double T_liquidus = T_solidus_liquidus(pressure, C_bar, false);

              std::vector<double> Tm = melting_temperatures(pressure);
              std::vector<double> K = partition_coefficients(pressure, std::max(T_solidus,std::min(T_liquidus,temperature)));
           
              double Pmax = 4.75e9;
              double Fmass_new = 0.0;
              double mcl = 0.0;
              double ccl = 0.0;
              double mcs = 0.0;
              double ccs = 0.0;
              if(this->get_adiabatic_conditions().pressure(in.position[q]) < Pmax)
              {
              // Calculate equilibrium melt fraction.
              //double Fvol_new = melt_fractions(Tm, K, C_bar, Fmass_old);
              //double Fvol_new = (rho_l/avg_rho)*Fmass_new; 
              //double Fmass_new = Fvol_new*(avg_rho/rho_l);

              // newton solver for calculating new mass fraction of melt.    
              Fmass_new = Fmass_old;
              int n      =  0;       
              const double r_tol = 1e-10;
              const int its_tol   = 100;
              double residual = 0.;
              for (unsigned int i=0; i<n_components; ++i)
                residual += C_bar[i]/(Fmass_old + (1-Fmass_old)*K[i]) - C_bar[i]/(Fmass_old/K[i] + (1-Fmass_old));

              if(temperature <= T_solidus)
                Fmass_new = 0;
              else if(temperature >= T_liquidus)
                Fmass_new = 1;
              else
              {
                while (abs(residual) > r_tol) 
                {
                  double dr_df = 0;
                  double term1 = 0;
                  double term2 = 0;
                  for (unsigned int i=0; i<n_components; ++i)
                  {
                    double numerator1 = C_bar[i] * (1.0 - K[i]);
                    double denominator1 = std::pow((Fmass_new + (1.0 - Fmass_new) * K[i]), 2);
                    term1 += numerator1 / denominator1;

                    double numerator2 = C_bar[i] * (1.0 / K[i] - 1.0);
                    double denominator2 = std::pow((Fmass_new / K[i] + (1.0 - Fmass_new)), 2);
                    term2 += numerator2 / denominator2;
                  }

                  dr_df = -term1+term2;

                  double a = 1;
                  while (((Fmass_new - a*residual/dr_df) < -1e-16) || (Fmass_new - a*residual/dr_df > 1-1e-16))
                  {
                    a = a/2;
                    if (a<1e-6)
                        AssertThrow(false, ExcMessage("a too small."));
                  }

                  Fmass_new = Fmass_new - a*residual/dr_df;

                  residual = 0.;
                  for (unsigned int i=0; i<n_components; ++i)
                    residual += C_bar[i]/(Fmass_new + (1-Fmass_new)*K[i]) - C_bar[i]/(Fmass_new/K[i] + (1-Fmass_new));

                  n = n+1;
                  if (n==its_tol)
                    AssertThrow(false, ExcMessage("No convergence"));
                }
              }

              // Calculate new Cl and Cs values, and limit all between 0 and 1.

              
              Fmass_new = std::max(0.0, std::min(1.0, Fmass_new));
              mcl = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new + (1 - Fmass_new) * K[1])));
              mcs = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new / K[1] + (1 - Fmass_new))));
              ccl = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new + (1 - Fmass_new) * K[2])));
              ccs = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new / K[2] + (1 - Fmass_new))));
              }
              else
              {
                  Fmass_new = 0.0;
                  mcl = 0.0;
                  ccl = 0.0;
                  mcs = 0.25;
                  ccs = 0.0005;
              }

              /*if (this->get_geometry_model().depth(in.position[q]) < 10e3)
              {
                  Fmass_new = 0.0;
                  mcl = 0.0;
                  ccl = 0.0;
                  mcs = 0.0;
                  ccs = 0.0;
              }*/

              // Calculate the reaction rates. I think this should be abs so its always positive?
              double Fvol_new = (rho_l/avg_rho)*Fmass_new;

              std::vector<double> Csf (n_components);
              std::vector<double> Clf (n_components);
              std::vector<double> CGamma (n_components);
              std::vector<double> Delta (n_components);
              std::vector<double> Gamma (n_components);
              std::vector<double> dcs (n_components);
              std::vector<double> dcl (n_components);

              // Calculate reaction rates ////////
              double R  =  rho_s/melting_time_scale;

              double GammaNet  =  R * (Fmass_new - Fmass_old);
 
              // Setup new compositions in order.
              std::vector<double> c_se = {(1 - mcs - ccs), mcs, ccs};
              std::vector<double> c_le = {(1 - mcl - ccl), mcl, ccl};
              for (unsigned int i=0; i<n_components; ++i)
              {
                Csf[i] = c_l[i]*K[i];
                Clf[i] = c_s[i]/K[i];

              if(GammaNet < 0)
                CGamma[i] = Csf[i];
              else if(GammaNet >= 0)
                CGamma[i] = Clf[i];

              Delta[i] = R*(Fmass_new*(c_le[i] - CGamma[i]) - Fmass_old*(c_l[i] - CGamma[i]));
              Gamma[i] = CGamma[i]*GammaNet + Delta[i];
              }

              double G2 = Gamma[0]+Gamma[1]+Gamma[2];

            for (unsigned int i=0; i<n_components; ++i)
            {
              dcs[i] = -(Gamma[i] - c_s[i]*G2)/(std::max(1e-6,(1 - Fmass_old))*rho_s);
              dcl[i] = (Gamma[i] - c_l[i]*G2)/(std::max(1e-6,(Fmass_old))*rho_s);
            }

              // Melt reaction rate using mass fraction
              double df = G2/rho_s;



              // WHAT TO OUTPUT: convert f back to volume fraction? It is mass fraction at the moment.
              // because depletion is a volume-based, and not a mass-based property that is advected,
              // additional scaling factors on the right hand side apply

              // We use the full compositions as the old values here.
              // so the rate includes any deviations outside of 0 and 1.
              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                  out.reaction_terms[q][c] = 0.0;
                  if (reaction_rate_out != nullptr && in.requests_property(MaterialProperties::reaction_rates) && this->get_timestep_number() > 0)
                  {
                    if (c == melt_idx)
                    {
                          // Melt reaction rate. Convert to volume for field. 
                          double rate = df*(rho_l/avg_rho); 
                          rate = std::max(rate*melting_time_scale, -in.composition[q][c]);
                          //double rate2 = std::max((Fvol_new-in.composition[q][c]), -in.composition[q][c]);

                          /*if (this->get_geometry_model().depth(in.position[q]) < 10e3)
                          {
                            Fvol_new = 0;
                            rate = std::max((Fvol_new-in.composition[q][c]), -in.composition[q][c]);
                          }*/

                          reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;

                    }
                    else if (c == mcs_idx)
                    {
                          // solid morb reaction rate.
                          double rate = std::max(dcs[1]*melting_time_scale, -in.composition[q][c]);
                          //double rate2 = std::max((mcs-in.composition[q][c]), -in.composition[q][c]);
                          reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                    }
                    else if (c == mcl_idx)
                    {
                          // liquid morb reaction rate.
                          double rate = std::max(dcl[1]*melting_time_scale, -in.composition[q][c]);
                          //double rate2 = std::max((mcl-in.composition[q][c]), -in.composition[q][c]);
                          reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                    }
                    else if (c == ccs_idx)
                    {
                          // solid cmorb reaction rate.
                          double rate = std::max(dcs[2]*melting_time_scale, -in.composition[q][c]);
                          //double rate2 = std::max((ccs-in.composition[q][c]), -in.composition[q][c]);
                          reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                    }
                    else if (c == ccl_idx)
                    {
                          // liquid cmorb reaction rate.
                          double rate = std::max(dcl[2]*melting_time_scale, -in.composition[q][c]);
                          //double rate2 = std::max((ccl-in.composition[q][c]), -in.composition[q][c]);
                          reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                    }
                    else
                      reaction_rate_out->reaction_rates[q][c] = 0.0;
                  }
                }

                out.entropy_derivative_pressure[q]    = 0.0;
                out.entropy_derivative_temperature[q] = 0.0;
            }
      }

template <int dim>
      void
      Co2Melt<dim>::
      calculate_fluid_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out,
                              const double reference_T) const
      {
        MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();

        if (melt_out != nullptr)
          {

            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                double porosity = std::max(0.0, std::min(in.composition[i][melt_idx],1.0));
                double morb_cl =  std::max(0.0, std::min(in.composition[i][mcl_idx],1.0));
                double cmorb_cl =  std::max(0.0, std::min(in.composition[i][ccl_idx],1.0));
                double dunite_cl = 1 - morb_cl - cmorb_cl;

                melt_out->fluid_viscosities[i] = viscosity_fluid*((pow(1.0,morb_cl))*(pow(10.0,dunite_cl))*(pow(0.01,cmorb_cl)));
                melt_out->permeabilities[i] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);

                // first, calculate temperature dependence of density
                double temperature_dependence = 1.0;
                if (this->include_adiabatic_heating ())
                  {
                    // temperature dependence is 1 - alpha * (T - T(adiabatic))
                    temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                              * out.thermal_expansion_coefficients[i];
                  }
                else
                  temperature_dependence -= (in.temperature[i] - reference_T) * out.thermal_expansion_coefficients[i];

                // the fluid compressibility includes two parts, a constant compressibility, and a pressure-dependent one
                // this is a simplified formulation, experimental data are often fit to the Birch-Murnaghan equation of state
                //const double fluid_compressibility = melt_compressibility / (1.0 + in.pressure[i] * melt_bulk_modulus_derivative * melt_compressibility);

                melt_out->fluid_densities[i] = rho_l * temperature_dependence;
                //reference_rho_fluid * std::exp(fluid_compressibility * (in.pressure[i] - this->get_surface_pressure()))
                                               //* temperature_dependence;

                melt_out->fluid_density_gradients[i] = 0.;

                const double phi_0 = 0.05;
                porosity = std::max(std::min(porosity,0.995),1e-4);
                melt_out->compaction_viscosities[i] = xi_0 * phi_0 / porosity;

                double visc_temperature_dependence = 1.0;
                if (this->include_adiabatic_heating ())
                  {
                    const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
                    visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
                  }
                else
                  {
                    const double delta_temp = in.temperature[i]-reference_T;
                    const double T_dependence = (thermal_bulk_viscosity_exponent == 0.0
                                                 ?
                                                 0.0
                                                 :
                                                 thermal_bulk_viscosity_exponent*delta_temp/reference_T);
                    visc_temperature_dependence = std::max(std::min(std::exp(-T_dependence),1e4),1e-4);
                  }
                melt_out->compaction_viscosities[i] *= visc_temperature_dependence;


              }
          }

        if (this->include_melt_transport() && in.requests_property(MaterialProperties::viscosity))
          {
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                const double porosity = std::min(1.0, std::max(in.composition[i][melt_idx],0.0));
                double viscosity  = out.viscosities[i]*std::exp(- alpha_phi * porosity);
                out.viscosities[i] = std::max(viscosity,5e15);
              }
          }
      }


      template <int dim>
      double
      Co2Melt<dim>::
      melt_fraction (std::vector<double> composition) const
      {
        // Here we return the current field.
        // TODO: should this call calculate_reaction_rates and return new field?
        double Fvol_old = std::max(0.0, std::min(composition[melt_idx],1.0));  
        return Fvol_old;
      }

      template <int dim>
      std::tuple<double, double, std::vector<double>, std::vector<double>>
      Co2Melt<dim>::
      equilibrium (std::vector<double> composition, 
                     const double temperature, 
                     const double pressure,
                     const double p2) const
      {
        // Define component dependent parameters 
        std::vector<double> C_bar (n_components);

        // Get the compositional values, and limit between 0 and 1.
        double morb_cl =  std::max(0.0, std::min(composition[mcl_idx],1.0));
        double morb_cs =  std::max(0.0, std::min(composition[mcs_idx],1.0));
        double cmorb_cl =  std::max(0.0, std::min(composition[ccl_idx],1.0));
        double cmorb_cs =  std::max(0.0, std::min(composition[ccs_idx],1.0));
        double Fvol_old =  std::max(0.0, std::min(composition[melt_idx],1.0));       

              // Calculate dunite and order liquid and solid components
              double dunite = 1 - morb_cs - cmorb_cs;
              double dunite_l = 1 - morb_cl - cmorb_cl;
              std::vector<double> c_s = {dunite, morb_cs, cmorb_cs};
              std::vector<double> c_l = {dunite_l, morb_cl, cmorb_cl};

              // We track the volume of melt, convert to mass here.
              double avg_rho = Fvol_old*rho_l + (1 - Fvol_old)*rho_s;
              double Fmass_old = Fvol_old*avg_rho/rho_l;

              // Now that things are ordered, find the bulk composition for each component.
              for (unsigned int i=0; i<n_components; ++i)
                C_bar[i] = Fmass_old*c_l[i] + (1-Fmass_old)*c_s[i];

              const double T_solidus = T_solidus_liquidus(pressure, C_bar, true);
              const double T_liquidus = T_solidus_liquidus(pressure, C_bar, false);

              std::vector<double> Tm = melting_temperatures(pressure);
              std::vector<double> K = partition_coefficients(pressure, std::max(T_solidus,std::min(T_liquidus,temperature)));
           
              double Pmax = 4.75e9;
              double Fmass_new = 0.0;
              double mcl = 0.0;
              double ccl = 0.0;
              double mcs = 0.0;
              double ccs = 0.0;
              if(p2 < Pmax)
              {
              // Calculate equilibrium melt fraction.
              //double Fvol_new = melt_fractions(Tm, K, C_bar, Fmass_old);
              //double Fvol_new = (rho_l/avg_rho)*Fmass_new; 
              //double Fmass_new = Fvol_new*(avg_rho/rho_l);

              // newton solver for calculating new mass fraction of melt.    
              Fmass_new = Fmass_old;
              int n      =  0;       
              const double r_tol = 1e-5;
              const int its_tol   = 100;
              double residual = 0.;
              for (unsigned int i=0; i<n_components; ++i)
                residual += C_bar[i]/(Fmass_old + (1-Fmass_old)*K[i]) - C_bar[i]/(Fmass_old/K[i] + (1-Fmass_old));

              if(temperature <= T_solidus)
                Fmass_new = 0;
              else if(temperature >= T_liquidus)
                Fmass_new = 1;
              else
              {
                while (abs(residual) > r_tol) 
                {
                  double dr_df = 0;
                  double term1 = 0;
                  double term2 = 0;
                  for (unsigned int i=0; i<n_components; ++i)
                  {
                    double numerator1 = C_bar[i] * (1.0 - K[i]);
                    double denominator1 = std::pow((Fmass_new + (1.0 - Fmass_new) * K[i]), 2);
                    term1 += numerator1 / denominator1;

                    double numerator2 = C_bar[i] * (1.0 / K[i] - 1.0);
                    double denominator2 = std::pow((Fmass_new / K[i] + (1.0 - Fmass_new)), 2);
                    term2 += numerator2 / denominator2;
                  }

                  dr_df = -term1+term2;

                  double a = 1;
                  while (((Fmass_new - a*residual/dr_df) < -1e-16) || (Fmass_new - a*residual/dr_df > 1-1e-16))
                  {
                    a = a/2;
                    if (a<1e-6)
                        AssertThrow(false, ExcMessage("a too small."));
                  }

                  Fmass_new = Fmass_new - a*residual/dr_df;

                  residual = 0.;
                  for (unsigned int i=0; i<n_components; ++i)
                    residual += C_bar[i]/(Fmass_new + (1-Fmass_new)*K[i]) - C_bar[i]/(Fmass_new/K[i] + (1-Fmass_new));

                  n = n+1;
                  if (n==its_tol)
                    AssertThrow(false, ExcMessage("No convergence"));
                }
              }

              // Calculate new Cl and Cs values, and limit all between 0 and 1.

              
              Fmass_new = std::max(0.0, std::min(1.0, Fmass_new));
              mcl = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new + (1 - Fmass_new) * K[1])));
              mcs = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new / K[1] + (1 - Fmass_new))));
              ccl = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new + (1 - Fmass_new) * K[2])));
              ccs = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new / K[2] + (1 - Fmass_new))));
              }
              else
              {
                  Fmass_new = 0.0;
                  mcl = 0.0;
                  ccl = 0.0;
                  mcs = 0.25;
                  ccs = 0.0005;
              }

              /*if (this->get_geometry_model().depth(in.position[q]) < 10e3)
              {
                  Fmass_new = 0.0;
                  mcl = 0.0;
                  ccl = 0.0;
                  mcs = 0.0;
                  ccs = 0.0;
              }*/

              // Calculate the reaction rates. I think this should be abs so its always positive?
              double Fvol_new = (rho_l/avg_rho)*Fmass_new;

              std::vector<double> Csf (n_components);
              std::vector<double> Clf (n_components);
              std::vector<double> CGamma (n_components);
              std::vector<double> Delta (n_components);
              std::vector<double> Gamma (n_components);
              std::vector<double> dcs (n_components);
              std::vector<double> dcl (n_components);

              // Calculate reaction rates ////////
              double R  =  rho_s/melting_time_scale;

              double GammaNet  =  R * (Fmass_new - Fmass_old);
 
              // Setup new compositions in order.
              std::vector<double> c_se = {(1 - mcs - ccs), mcs, ccs};
              std::vector<double> c_le = {(1 - mcl - ccl), mcl, ccl};
              for (unsigned int i=0; i<n_components; ++i)
              {
                Csf[i] = c_l[i]*K[i];
                Clf[i] = c_s[i]/K[i];

              if(GammaNet < 0)
                CGamma[i] = Csf[i];
              else if(GammaNet >= 0)
                CGamma[i] = Clf[i];

              Delta[i] = R*(Fmass_new*(c_le[i] - CGamma[i]) - Fmass_old*(c_l[i] - CGamma[i]));
              Gamma[i] = CGamma[i]*GammaNet + Delta[i];
              }

              double G2 = Gamma[0]+Gamma[1]+Gamma[2];

            for (unsigned int i=0; i<n_components; ++i)
            {
              dcs[i] = -(Gamma[i] - c_s[i]*G2)/(std::max(1e-6,(1 - Fmass_old))*rho_s);
              dcl[i] = (Gamma[i] - c_l[i]*G2)/(std::max(1e-6,(Fmass_old))*rho_s);
            }

              // Melt reaction rate using mass fraction
              double df = G2/rho_s;                                                       


      // Convert melt parameters to volume.
      //return {Fmass_new*(rho_l/avg_rho), df*(rho_l/avg_rho), dcs, dcl};
      return {Fmass_new*(rho_l/avg_rho), df, dcs, dcl};
      }

      template <int dim>
      double
      Co2Melt<dim>::
      T_solidus_liquidus (const double pressure, 
                          std::vector<double> composition, 
                          bool compute_solidus) const
      {
        // TODO: Exclude invalid compositions (that do not sum up to 1)?
        const std::vector<double> Tm = melting_temperatures(pressure);

        // Set starting guess for Tsol
        const double minTm = *std::min_element(Tm.begin(), Tm.end());
        const double maxTm = *std::max_element(Tm.begin(), Tm.end());
        double mean_Tm = 0.0;

        for (unsigned int i=0; i<n_components; ++i) 
        {
            double comp = composition[i];
            if (n_components == 0)
                comp = 1 - composition[i+1]; 

            mean_Tm += comp * Tm[i];
        }

        double T_solidus = std::max(minTm, std::min(maxTm, mean_Tm));

        std::vector<double> K = partition_coefficients(pressure, T_solidus);

        // Get residual for sum(ci_bar/Ki) = 1 or sum(ci_bar*Ki) = 1 (Equations 8 + 9)
        double residual = compute_residual(composition, K, compute_solidus);

        unsigned int n                    =  0;     // initialize iteration count
        const double tolerance            =  1e-10; //1e-10; // tolerance for Newton residual
        const unsigned int max_iterations =  2000;   // maximum number of iterations
        const double eps_T                =  5;     // temperature perturbation for finite differencing, degrees

        while (std::abs(residual) > tolerance) 
        {
          // Compute partition coefficients Ki at T+eps_T
          K = partition_coefficients(pressure, T_solidus + eps_T);

          // Get residual at T + eps_T
          double residual_plus_eps_T = compute_residual(composition, K, compute_solidus);

          // Compute partition coefficients Ki at T-eps_T
          K = partition_coefficients(pressure, T_solidus - eps_T);

          // Get residual at T + eps_T
          double residual_minus_eps_T = compute_residual(composition, K, compute_solidus);

          // Finite difference drdT = (r(T+eps_T)-r(T-eps_T))/2/eps_T
          const double dresidualdT  =  (residual_plus_eps_T - residual_minus_eps_T) / (2 * eps_T);

          // Apply Newton correction to current guess of Tsol
          T_solidus = T_solidus - 0.5 * residual/dresidualdT;

          // Compute partition coefficients Ki at Tsol
          K = partition_coefficients(pressure, T_solidus);

          // Get residual at T_solidus
          residual = compute_residual(composition, K, compute_solidus);

          ++n;

          if (n == max_iterations) 
          {
            std::cerr << "!!! Newton solver for solidus/liquidus T has not converged after " << residual << " iterations !!!" << std::endl;
            break;
          }
        }
      return T_solidus;
    }

    template <int dim>
    std::vector<double>
    Co2Melt<dim>::melting_temperatures(const double pressure) const
    {
        std::vector<double> Tm (n_components);

        const double Pmax = 6e9;
        if (pressure <= Pmax)
          for (unsigned int i=0; i<n_components; ++i) 
            Tm[i]  =  T0[i] + A[i] * pressure + B[i] * pressure * pressure;

        else
        {
          // safeguard: continue melting point with linear slope above Pmax
          const double dP = 1e7;
          for (unsigned int i=0; i<n_components; ++i) 
          {
            bool ind = pressure > Pmax;
            const double T0_at_Pmax = T0[i] + A[i] * Pmax + B[i] * Pmax * Pmax;
            const double dTdP = ((A[i]*Pmax + B[i]*Pmax*Pmax) - (A[i]*(Pmax-1e7) + B[i]*(Pmax-1e7)*(Pmax-1e7)))/dP;
            Tm[i] = T0_at_Pmax + dTdP * (pressure-Pmax);
          }
        }
        return Tm;
    }


    template <int dim>
    std::vector<double>
    Co2Melt<dim>::partition_coefficients(const double pressure, 
                                                const double temperature) const
    {
        std::vector<double> K (n_components);

        std::vector<double> Tm = melting_temperatures(pressure);

        // Parameterization after Rudge, Bercovici, & Spiegelman (2010)
        for (unsigned int i=0; i<n_components; ++i) 
        {
          double Ls = L[i]/T0[i]*temperature;
            K[i] = std::exp(Ls/R[i] * (1./temperature - 1./Tm[i]));
        }

        return K;
    }

    template <int dim>
    double 
    Co2Melt<dim>::compute_residual (std::vector<double> composition,
                      std::vector<double> K,
                      bool compute_solidus) const
    {
      double residual = -1.;
      for (unsigned int i=0; i<n_components; ++i) 
        if (compute_solidus)
          residual += composition[i]/K[i]; // solidus
        else
          residual += composition[i]*K[i]; // liquidus

      return residual;
    }     


      template <int dim>
      void
      Co2Melt<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Reaction model");
          {
          prm.enter_subsection("Co2 model");
          {
              prm.declare_entry ("T0", "1085.7",
                                 Patterns::List(Patterns::Double (0.)),
                                 "Constant parameter in the quadratic "
                                 "function that approximates the solidus "
                                 "of peridotite. "
                                 "Units: \\si{\\degreeCelsius}.");
              prm.declare_entry ("A", "1.329e-7",
                                 Patterns::List(Patterns::Double (0.)),
                                 "Prefactor of the linear pressure term "
                                 "in the quadratic function that approximates "
                                 "the solidus of peridotite. "
                                 "\\si{\\degreeCelsius\\per\\pascal}.");
              prm.declare_entry ("B", "-5.1e-18",
                                 Patterns::List(Patterns::Double ()),
                                 "Prefactor of the quadratic pressure term "
                                 "in the quadratic function that approximates "
                                 "the solidus of peridotite. "
                                 "\\si{\\degreeCelsius\\per\\pascal\\squared}.");
              prm.declare_entry ("L", "1475.0",
                                 Patterns::List(Patterns::Double (0.)),
                                 "Constant parameter in the quadratic "
                                 "function that approximates the lherzolite "
                                 "liquidus used for calculating the fraction "
                                 "of peridotite-derived melt. "
                                 "Units: \\si{\\degreeCelsius}.");
              prm.declare_entry ("R", "8.0e-8",
                                 Patterns::List(Patterns::Double (0.)),
                                 "Prefactor of the linear pressure term "
                                 "in the quadratic function that approximates "
                                 "the  lherzolite liquidus used for "
                                 "calculating the fraction of peridotite-"
                                 "derived melt. "
                                 "\\si{\\degreeCelsius\\per\\pascal}.");
              prm.declare_entry ("Number of components", "3",
                                 Patterns::Integer(),
                                 "Prefactor of the linear pressure term "
                                 "in the quadratic function that approximates "
                                 "the  lherzolite liquidus used for "
                                 "calculating the fraction of peridotite-"
                                 "derived melt. "
                                 "\\si{\\degreeCelsius\\per\\pascal}.");            

              prm.declare_entry ("Fluid density", "2700",
                    Patterns::List(Patterns::Double (0.)),
                    "Constant parameter in the quadratic "
                    "function that approximates the solidus "
                    "of peridotite. "
                    "Units: \\si{\\degreeCelsius}.");

              prm.declare_entry ("Solid density", "3200",
                    Patterns::List(Patterns::Double (0.)),
                    "Constant parameter in the quadratic "
                    "function that approximates the solidus "
                    "of peridotite. "
                    "Units: \\si{\\degreeCelsius}.");
              prm.declare_entry ("Melting time scale for operator splitting", "1e2",
                    Patterns::Double (0.),
                    "Because the operator splitting scheme is used, the porosity field can not "
                    "be set to a new equilibrium melt fraction instantly, but the model has to "
                    "provide a melting time scale instead. This time scale defines how fast melting "
                    "happens, or more specifically, the parameter defines the time after which "
                    "the deviation of the porosity from the equilibrium melt fraction will be "
                    "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                    "to the time step size, the reaction will be so fast that the porosity is very "
                    "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                    "if the melting time scale is large compared to the time step size, almost no "
                    "melting and freezing will occur."
                    "\n\n"
                    "Also note that the melting time scale has to be larger than or equal to the reaction "
                    "time step used in the operator splitting scheme, otherwise reactions can not be "
                    "computed. "
                    "Units: yr or s, depending on the ``Use years in output instead of seconds'' parameter.");
              prm.declare_entry ("Reference bulk viscosity", "1e22",
                                Patterns::Double (0.),
                                "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                                "This viscosity may be modified by both temperature and porosity "
                                "dependencies. Units: \\si{\\pascal\\second}.");
              prm.declare_entry ("Reference melt viscosity", "10.",
                                Patterns::Double (0.),
                                "The value of the constant melt viscosity $\\viscosity_fluid$. Units: \\si{\\pascal\\second}.");
              prm.declare_entry ("Exponential melt weakening factor", "27.",
                                Patterns::Double (0.),
                                "The porosity dependence of the viscosity. Units: dimensionless.");
              prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                                Patterns::Double (0.),
                                "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                                "See the general documentation "
                                "of this model for a formula that states the dependence of the "
                                "viscosity on this factor, which is called $\\beta$ there.");
              prm.declare_entry ("Melt compressibility", "0.0",
                                Patterns::Double (0.),
                                "The value of the compressibility of the melt. "
                                "Units: \\si{\\per\\pascal}.");
              prm.declare_entry ("Melt bulk modulus derivative", "0.0",
                                Patterns::Double (0.),
                                "The value of the pressure derivative of the melt bulk "
                                "modulus. "
                                "Units: None.");
              prm.declare_entry ("Reference permeability", "1e-6",
                                Patterns::Double(),
                                "Reference permeability of the solid host rock."
                                "Units: \\si{\\meter\\squared}.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
      
      


      template <int dim>
      void
      Co2Melt<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Reaction model");
        {
        prm.enter_subsection("Co2 model");
        {
            n_components = prm.get_integer("Number of components");
            T0 = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("T0"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            A = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("A"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");                                                              
            B = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("B"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            L = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("L"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            R = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("R"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            rho_l         = prm.get_double ("Fluid density");
            rho_s         = prm.get_double ("Solid density");
            melting_time_scale         = prm.get_double ("Melting time scale for operator splitting");

            if (this->convert_output_to_years() == true)
              melting_time_scale *= year_in_seconds;

            // Get ID for all fields. Should I do this once here or do it in the main function?
            AssertThrow(this->introspection().compositional_name_exists("porosity"), ExcMessage("A porosity field is needed to use the co2 plugin."));
            melt_idx = this->introspection().compositional_index_for_name("porosity");

            AssertThrow(this->introspection().compositional_name_exists("morb_cl"), ExcMessage("A morb_cl field is needed to use the co2 plugin."));
            mcl_idx = this->introspection().compositional_index_for_name("morb_cl");

            AssertThrow(this->introspection().compositional_name_exists("morb_cs"), ExcMessage("A morb_cs field is needed to use the co2 plugin."));
            mcs_idx = this->introspection().compositional_index_for_name("morb_cs");
          
            AssertThrow(this->introspection().compositional_name_exists("cmorb_cl"), ExcMessage("A cmorb_cl field is needed to use the co2 plugin."));
            ccl_idx = this->introspection().compositional_index_for_name("cmorb_cl");

            AssertThrow(this->introspection().compositional_name_exists("cmorb_cs"), ExcMessage("A cmorb_cs field is needed to use the co2 plugin."));
            ccs_idx = this->introspection().compositional_index_for_name("cmorb_cs");

            xi_0                       = prm.get_double ("Reference bulk viscosity");
            viscosity_fluid            = prm.get_double ("Reference melt viscosity");
            thermal_bulk_viscosity_exponent = prm.get_double ("Thermal bulk viscosity exponent");
            alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
            melt_compressibility       = prm.get_double ("Melt compressibility");
            melt_bulk_modulus_derivative = prm.get_double ("Melt bulk modulus derivative");
            reference_permeability     = prm.get_double ("Reference permeability");
            }
            prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace ReactionModel \
  { \
    template class Co2Melt<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
