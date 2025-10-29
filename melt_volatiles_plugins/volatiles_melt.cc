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

#include </mnt/vast-nhr/home/derekjohn.neuharth/u16318/software/co2/aspect/melt_volatiles_plugins/volatiles_melt.h>
#include <aspect/utilities.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <deal.II/base/parameter_handler.h>
#include <chrono>

#include <aspect/simulator.h>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {


      template <int dim>
      VolatilesMelt<dim>::VolatilesMelt()
        = default;

      template <int dim>
      double
      VolatilesMelt<dim>::
      reference_darcy_coefficient () const
      {
        // 0.01 = 1% melt
        // Darcy coefficient should vary with melt viscosity.
        return reference_permeability * Utilities::fixed_power<3>(0.01) / viscosity_fluid;
      }

      template <int dim>
      void
      VolatilesMelt<dim>::
      calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                      typename Interface<dim>::MaterialModelOutputs &out) const
      {
          // Reaction rates needed for operator splitting. This model doesn't consider a case
          // where there is no operator splitting that uses reaction terms instead.
          const std::shared_ptr<ReactionRateOutputs<dim>>
          reaction_rate_out = out.template get_additional_output_object<ReactionRateOutputs<dim>>();

          // Enthalpy outputs for the latent heat mel plugin.
          const std::shared_ptr<EnthalpyOutputs<dim>> enthalpy_out 
            = out.template get_additional_output_object<EnthalpyOutputs<dim>>();

          double reaction_time_step_size = 1.0;
          if (this->simulator_is_past_initialization())
          {
            const unsigned int number_of_reaction_steps = std::max(static_cast<unsigned int>(this->get_timestep() / this->get_parameters().reaction_time_step),
                                                                  std::max(this->get_parameters().reaction_steps_per_advection_step,1U));
            reaction_time_step_size = this->get_timestep() / static_cast<double>(number_of_reaction_steps);
          }

          for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
            {
              const double temperature = in.temperature[q];
              const double pressure = this->get_adiabatic_conditions().pressure(in.position[q]) > 101325.
              ? 
              this->get_adiabatic_conditions().pressure(in.position[q])
              :
              101325.;

              std::vector<double> composition(this->n_compositional_fields());
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                      composition[c] = in.composition[q][c];

              const double ycord = in.position[q](1);
              const double xcord = in.position[q](0);
              const double depth = this->get_geometry_model().depth(in.position[q]);
              double porosity = in.composition[q][porosity_idx];
              // Ignore melt fraction, get melt reaction rate (volume)
              // and solid and liquid reaction rates, ordered as dunite (background field), morb, cmorb, hmorb.
              // Note: At the moment compositions are hardcoded in assuming there is always 4 components.
              const double rho_s = out.densities[q]; //3200
              auto [vfrac, melt_reaction_rate, solid_reaction_rates, liquid_reaction_rates, enthalpy] = equilibrium(composition, temperature, pressure, ycord, rho_s, q, xcord);

              for (unsigned int c=0; c<in.composition[q].size(); ++c)
                {
                
                  if (reaction_rate_out != nullptr && in.requests_property(MaterialProperties::reaction_rates))
                  {
                    if (c == porosity_idx)
                    {
                          reaction_rate_out->reaction_rates[q][c] = get_reaction_rate(in.composition[q][c], 
                                                                                      melt_reaction_rate, 
                                                                                      reaction_time_step_size,
                                                                                      depth,
                                                                                      xcord);
                    }
                    else if (c == mcs_idx)
                    {
                          // solid morb reaction rate, between zero and 1. Because the components 
                          reaction_rate_out->reaction_rates[q][c] = get_reaction_rate(in.composition[q][c], 
                                                          solid_reaction_rates[1], 
                                                          reaction_time_step_size,
                                                          depth,
                                                          xcord);

                    }
                    else if (c == mcl_idx)
                    {
                          // liquid morb reaction rate.
                          reaction_rate_out->reaction_rates[q][c] = get_reaction_rate(in.composition[q][c], 
                                                          liquid_reaction_rates[1], 
                                                          reaction_time_step_size,
                                                          depth,
                                                          xcord);
                    }
                    else if (c == ccs_idx)
                    {
                          // solid cmorb reaction rate.
                          reaction_rate_out->reaction_rates[q][c] = get_reaction_rate(in.composition[q][c], 
                                                          solid_reaction_rates[2], 
                                                          reaction_time_step_size,
                                                          depth,
                                                          xcord);                            
                    }
                    else if (c == ccl_idx)
                    {
                          // liquid cmorb reaction rate.
                          reaction_rate_out->reaction_rates[q][c] = get_reaction_rate(in.composition[q][c], 
                                                          liquid_reaction_rates[2], 
                                                          reaction_time_step_size,
                                                          depth,
                                                          xcord);
                    }
                    else if (c == hcs_idx)
                    {
                          // solid cmorb reaction rate.
                          reaction_rate_out->reaction_rates[q][c] = get_reaction_rate(in.composition[q][c], 
                                                          solid_reaction_rates[3], 
                                                          reaction_time_step_size,
                                                          depth,
                                                          xcord);
                    }
                    else if (c == hcl_idx)
                    {
                      // liquid cmorb reaction rate.
                      reaction_rate_out->reaction_rates[q][c] = get_reaction_rate(in.composition[q][c], 
                                                      liquid_reaction_rates[3], 
                                                      reaction_time_step_size,
                                                      depth,
                                                      xcord);
                    }
                    else
                      reaction_rate_out->reaction_rates[q][c] = 0.0;
                  }
                }

                if (enthalpy_out != nullptr)
                  enthalpy_out->enthalpies_of_fusion[q] = enthalpy;    

                out.entropy_derivative_pressure[q]    = 0.0;
                out.entropy_derivative_temperature[q] = 0.0;
            }
      }

template <int dim>
      void
      VolatilesMelt<dim>::
      calculate_fluid_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                              typename Interface<dim>::MaterialModelOutputs &out,
                              const double reference_T) const
      {
        MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();

        // Next, find fluid outputs.
        if (melt_out != nullptr)
          {
            for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
              {
                double porosity = std::max(0.0, std::min(in.composition[i][porosity_idx],1.0));
                double morb_cl =  std::max(0.0, std::min(in.composition[i][mcl_idx],1.0));
                double cmorb_cl =  std::max(0.0, std::min(in.composition[i][ccl_idx],1.0));
                double hmorb_cl =  std::max(0.0, std::min(in.composition[i][hcl_idx],1.0));

                // Near the surface we sometimes get liquid values above one.
                // here we force them to 1 and assume there is no dunite.
                double comp_sum = morb_cl + cmorb_cl + hmorb_cl;
                if(comp_sum > 1)
                {
                  morb_cl = morb_cl/comp_sum;
                  cmorb_cl = cmorb_cl/comp_sum;
                  hmorb_cl = hmorb_cl/comp_sum;
                }

                double dunite_cl = std::max(0. ,(1 - morb_cl - cmorb_cl - hmorb_cl));
                
                // We should maybe double check that the liquid components sum to unity.
                melt_out->fluid_viscosities[i] = viscosity_fluid*((pow(1.0, morb_cl))*(pow(10.0, dunite_cl))*(pow(0.1, hmorb_cl))*(pow(0.01, cmorb_cl)));
            
                melt_out->permeabilities[i] = std::min(reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity), maximum_permeability);

                // At the moment we don't include compressibility.
                melt_out->fluid_densities[i] = out.densities[i] - fluid_density_difference;
                melt_out->fluid_density_gradients[i] = 0.;

                // limit porosity to disaggregation threshold
                porosity = std::min(0.3, porosity);
                const double porosity_threshold = 0.01 * std::pow(this->get_melt_handler().melt_parameters.melt_scaling_factor_threshold, 1./3.);

                // At the moment we use a defined reference solid viscosity. At some point should we use the model viscosity?
                // In this case, it probably shouldn't include changes related to the porosity. Move viscosity calculation below this?
                double viscosity = xi_0;
                if(this->get_timestep_number() > 0)
                  viscosity = out.viscosities[i];

                melt_out->compaction_viscosities[i] = compaction_viscosity_ratio * viscosity / std::max(porosity, porosity_threshold);
              }
          }

        // Find new viscosity.
        if (this->include_melt_transport() && in.requests_property(MaterialProperties::viscosity))
        {
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
            // cutoff for viscosity at 30%
            const double porosity = std::min(0.3, std::max(in.composition[i][porosity_idx],0.0));
            out.viscosities[i] = std::max(out.viscosities[i] * std::exp(- alpha_phi * porosity),1e15);
            }
        }
      }


      template <int dim>
      std::tuple<double, double, std::vector<double>, std::vector<double>, double>
      VolatilesMelt<dim>::
      equilibrium (std::vector<double> composition, 
                     const double temperature, 
                     const double pressure,
                     const double ycord,
                     const double rho_s,
                     const int ep,
                     const double x) const
      {
        // Define component dependent parameters 
        std::vector<double> C_bar (n_components);

        // Get the compositional values, and limit between 0 and 1.
        double morb_cl =  std::max(0.0, std::min(composition[mcl_idx],1.0));
        double morb_cs =  std::max(0.0, std::min(composition[mcs_idx],1.0));
        double cmorb_cl =  std::max(0.0, std::min(composition[ccl_idx],1.0));
        double cmorb_cs =  std::max(0.0, std::min(composition[ccs_idx],1.0));
        double hmorb_cl =  std::max(0.0, std::min(composition[hcl_idx],1.0));
        double hmorb_cs =  std::max(0.0, std::min(composition[hcs_idx],1.0));
        double Fvol_old =  std::max(0.0, std::min(composition[porosity_idx],1.0));
        const double rho_l = rho_s - fluid_density_difference;

        double reaction_time_step_size = 1.0;
        if (this->simulator_is_past_initialization())
        {
          const unsigned int number_of_reaction_steps = std::max(static_cast<unsigned int>(this->get_timestep() / this->get_parameters().reaction_time_step),
                                                                 std::max(this->get_parameters().reaction_steps_per_advection_step,1U));
          reaction_time_step_size = this->get_timestep() / static_cast<double>(number_of_reaction_steps);
        }

        // Near the surface we sometimes get liquid values above one.
        // here we force them to 1 and assume there is no dunite
        // Would there be a better way to do this?
        double comp_sum = morb_cl + cmorb_cl + hmorb_cl;
        if(comp_sum > 1)
        {
          morb_cl = morb_cl/comp_sum;
          cmorb_cl = cmorb_cl/comp_sum;
          hmorb_cl = hmorb_cl/comp_sum;
        }

        comp_sum = morb_cs + cmorb_cs + hmorb_cs;
        if(comp_sum > 1)
        {
          morb_cs = morb_cs/comp_sum;
          cmorb_cs = cmorb_cs/comp_sum;
          hmorb_cs = hmorb_cs/comp_sum;
        }

        // Find dunite and make sure it isn't below zero.
        double dunite = std::max(0. ,(1 - morb_cs - cmorb_cs - hmorb_cs)); 
        double dunite_l = std::max(0. ,(1 - morb_cl - cmorb_cl - hmorb_cl));    
          
        std::vector<double> c_s = {dunite, morb_cs, cmorb_cs, hmorb_cs};
        std::vector<double> c_l = {dunite_l, morb_cl, cmorb_cl, hmorb_cl};

        double avg_rho = Fvol_old*rho_l + (1. - Fvol_old)*rho_s;
        double avg_rho_new = avg_rho; // Will be updated later.
        double Fmass_old = Fvol_old*rho_l/avg_rho;

        if(ep == 1)
        {
          double cppm = (Fmass_old * cmorb_cl + (1 - Fmass_old)*cmorb_cs) * 20/100 * 1e6;
          double cmass = (Fmass_old * cmorb_cl * rho_l + (1 - Fmass_old)*cmorb_cs*rho_s) * 20/100;
          std::cout<<"Old: "<<cppm<<" "<<cmass<<std::endl;
        }
   
        // Now that things are ordered, find the bulk composition for each component.
        for (unsigned int i=0; i<n_components; ++i)
            C_bar[i] = Fmass_old*c_l[i] + (1-Fmass_old)*c_s[i];
      
        // Define parameters that will be returned.
        double Fmass_new = 0.0;
        std::vector<double> solid_reaction_rates (n_components, 0.0);
        std::vector<double> liquid_reaction_rates (n_components, 0.0);
        double melt_reaction_rate = 0.0;
        double enthalpy = 0.0;

        std::vector<double> T0 (n_components);
        std::vector<double> A (n_components);
        std::vector<double> B (n_components);
        std::vector<double> L (n_components);
        std::vector<double> R (n_components);

        initialize_component_values(A, B, L, T0, R, pressure);

        // Calculate the equilibrium values and reaction rates
        // if we are below the maximum solidus pressure.
        if(pressure < pressure_max)
        {
          const double T_solidus = T_solidus_liquidus(pressure, C_bar, true, A, B, L, T0, R);
          const double T_liquidus = T_solidus_liquidus(pressure, C_bar, false, A, B, L, T0, R);

          std::vector<double> Tm = melting_temperatures(pressure, A, B, L, T0);
          std::vector<double> K = partition_coefficients(pressure, std::max(T_solidus,std::min(T_liquidus,temperature)), A, B, L, T0, R);
          
          // Calculate equilibrium melt fraction.
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

          // Provide maximum limit to porosity.
          //if(Fmass_new*(avg_rho/rho_l) > 0.3)
          //  Fmass_new = 0.3 * (rho_l / avg_rho);

          double dcl = std::max(0.0, std::min(1.0, C_bar[0] / (Fmass_new + (1 - Fmass_new) * K[0])));
          double mcl = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new + (1 - Fmass_new) * K[1])));
          double ccl = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new + (1 - Fmass_new) * K[2])));     
          double hcl = std::max(0.0, std::min(1.0, C_bar[3] / (Fmass_new + (1 - Fmass_new) * K[3])));
          
          // Solid values, these aren't actually used for the reaction rates so can remove.
          double dcs = std::max(0.0, std::min(1.0, C_bar[0] / (Fmass_new / K[0] + (1 - Fmass_new))));
          double mcs = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new / K[1] + (1 - Fmass_new))));
          double ccs = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new / K[2] + (1 - Fmass_new))));
          double hcs = std::max(0.0, std::min(1.0, C_bar[3] / (Fmass_new / K[3] + (1 - Fmass_new))));
        
          // Define reaction rate parameters and calculate rates for each component.
          std::vector<double> Gamma (n_components);
          double GammaSum = 0.0;

          // In the paper they use the constant reference density, resulting in a 
          // constant R factor of 3. We use the model density, so there may be
          // some variation. How important is this? Maybe we don't want to use
          // model density as it will take into consideration other compositions.
          double R  =  avg_rho/melting_time_scale;

          // Check unity, setup new equilibirum liquid in order, and calculate reaction rates.
          comp_sum = dcl + ccl + mcl + hcl;
          dcl = dcl/comp_sum;
          mcl = mcl/comp_sum;
          ccl = ccl/comp_sum;
          hcl = hcl/comp_sum;

          // Check unity, setup new equilibirum liquid in order, and calculate reaction rates.
          comp_sum = dcs + ccs + mcs + hcs;
          dcs = dcs/comp_sum;
          mcs = mcs/comp_sum;
          ccs = ccs/comp_sum;
          hcs = hcs/comp_sum;

          std::vector<double> c_leq = {dcl, mcl, ccl, hcl};
          std::vector<double> c_seq = {dcs, mcs, ccs, hcs};
          if(use_fractional_melting)
          {
            std::vector<double> Csf (n_components);
            std::vector<double> Clf (n_components);
            std::vector<double> CGamma (n_components);
            std::vector<double> Delta (n_components);
            double GammaNet  =  R * (Fmass_new - Fmass_old);
            double gsum = 0;
            for (unsigned int i=0; i<n_components; ++i)
            {
              Csf[i] = c_l[i]*K[i];
              Clf[i] = c_s[i]/K[i];

              if(GammaNet < 0)
                CGamma[i] = Csf[i];
              else if(GammaNet >= 0)
                CGamma[i] = Clf[i];
            }

            // Sum CGamma's to make unity continuing.
            gsum = CGamma[0]+CGamma[1]+CGamma[2]+CGamma[3];
            for (unsigned int i=0; i<n_components; ++i)
            {
              CGamma[i] = CGamma[i]/gsum;

              // In paper, they mention you can use a different melt timescale for R.
              Delta[i] = R*(Fmass_new*(c_leq[i] - CGamma[i]) - Fmass_old*(c_l[i] - CGamma[i]));
              Gamma[i] = CGamma[i]*GammaNet + Delta[i];

              double Ls = L[i]/T0[i]*temperature;
              enthalpy += (Gamma[i]*Ls);

              // In matlab code, they also include a pressure and temperature change, must these be added
              // somewhere? One componenet of temperature change is in the latent heat/enthalpy and is included,
              // others I am not sure about.
              GammaSum += Gamma[i];
            }
          }
          else // Using batch melting.
          {
            for (unsigned int i=0; i<n_components; ++i)
            {
              Gamma[i] = R*(Fmass_new*c_leq[i] - Fmass_old*c_l[i]);

              double Ls = L[i]/T0[i]*temperature;
              enthalpy += (Gamma[i]*Ls);

              GammaSum += Gamma[i];
            }
          }

          // Now that we have GammaSum, find the solid and liquid reaction rates.
          // From eq. 17b and 17c in Keller and Katz, 2016
          for (unsigned int i=0; i<n_components; ++i)
          {
            solid_reaction_rates[i] = -(Gamma[i] - c_s[i]*GammaSum) / (std::max(1e-6,(1 - Fmass_new))*avg_rho);
            liquid_reaction_rates[i] = (Gamma[i] - c_l[i]*GammaSum) / (std::max(1e-6,Fmass_new)*avg_rho);
          }

          // Melt reaction rate using mass fraction
          melt_reaction_rate = GammaSum/avg_rho;

          // Find the mass fraction of melt we will have at the end of the reaction step. 
          double melt_reaction_step = Fmass_old + melt_reaction_rate * reaction_time_step_size;

          // Find the updated average density.
          avg_rho_new = rho_s / (1 - melt_reaction_step * (1 - rho_s / rho_l));

           // Here we find what melt value we will reach by the end of the reaction timestep,
           // and adjust the liquid component so the bulk composition is convserved. Only do this
           // if we have some melt.
          if(Fmass_new > 0)
          {
            for (unsigned int i=0; i<n_components; ++i)
            {
            // Find the equilibrium bulk composition, defined as cl + cs at the end of rhe reaction time step.
            double cl_reaction_step = melt_reaction_step * (c_l[i] + liquid_reaction_rates[i] * reaction_time_step_size);
            double cs_reaction_step = (1 - melt_reaction_step) * (c_s[i] + solid_reaction_rates[i] * reaction_time_step_size);
            double Cbar_reaction_step = cl_reaction_step + cs_reaction_step;

            // Find the difference in bulk composition between the equilibrium
            // value and where we will be at the end of the reaction step.
            double bulk_composition_change = Cbar_reaction_step - C_bar[i];

            // Adjust the liquid component to conserve volatiles.
            double liquid_rate_correction = 0;
            if(melt_reaction_step > 0)
              liquid_rate_correction = (bulk_composition_change)/(reaction_time_step_size*melt_reaction_step);

            // Change reaction rate so it is between 0 and 1.
            double update_liquid = liquid_reaction_rates[i]-liquid_rate_correction;
            if(c_l[i]+update_liquid*reaction_time_step_size < 0)
              liquid_reaction_rates[i] = -c_l[i]/reaction_time_step_size;
            else if(c_l[i]+update_liquid*reaction_time_step_size > 1)
              liquid_reaction_rates[i] = (1.0-c_l[i])/reaction_time_step_size;
            else
              liquid_reaction_rates[i] = update_liquid;
            }
          }

        if(ep == 1)
        {
          double cppm = (Fmass_new * ccl + (1 - Fmass_new)*ccs) * 20/100 * 1e6;
          double cmass = (Fmass_new * ccl * rho_l + (1 - Fmass_new)*ccs*rho_s) * 20/100;
          std::cout<<"New: "<<Fmass_old<<" "<<Fmass_new<<" "<<cppm<<" "<<cmass<<std::endl;
        }

          if(reaction_time_step_size > 0)
            melt_reaction_rate += Fvol_old * rho_l * (1.0 / avg_rho - 1.0 / avg_rho_new) / reaction_time_step_size;


          if(GammaSum != 0)
            enthalpy = enthalpy/GammaSum;         
      }   

      // Return values, with melt_fractions converted from mass fraction to volume fraction.
      return {Fmass_new*(avg_rho_new/rho_l), melt_reaction_rate*(avg_rho_new/rho_l), solid_reaction_rates, liquid_reaction_rates, enthalpy};
      }

      template <int dim>
      double
      VolatilesMelt<dim>::
      T_solidus_liquidus (const double pressure, 
                          std::vector<double> composition, 
                          bool compute_solidus,
                          std::vector<double> A,
                          std::vector<double> B,
                          std::vector<double> L,
                          std::vector<double> T0,
                          std::vector<double> R) const
      {
        // TODO: Exclude invalid compositions (that do not sum up to 1)?
        const std::vector<double> Tm = melting_temperatures(pressure, A, B, L, T0);

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

        std::vector<double> K = partition_coefficients(pressure, T_solidus, A, B, L, T0, R);

const double tolerance = 1e-10;
const unsigned int max_iterations = 200;
double residual = compute_residual(composition, K, compute_solidus);

//auto start = std::chrono::high_resolution_clock::now();
unsigned int n = 0;
while (std::abs(residual) > tolerance && n < max_iterations)
{
    // Adaptive perturbation
    const double eps_T = std::max(1e-6, 1e-6 * std::abs(T_solidus));

    // Central difference derivative
    K = partition_coefficients(pressure, T_solidus + eps_T, A, B, L, T0, R);
    double residual_plus = compute_residual(composition, K, compute_solidus);

    K = partition_coefficients(pressure, T_solidus - eps_T, A, B, L, T0, R);
    double residual_minus = compute_residual(composition, K, compute_solidus);

    const double dresidualdT = (residual_plus - residual_minus) / (2.0 * eps_T);

    if (std::abs(dresidualdT) < 1e-14) {
        std::cerr << "Derivative vanished, aborting Newton at iteration " << n << std::endl;
        break;
    }

    // Newton step with optional damping
    double delta_T = -residual / dresidualdT;
    T_solidus += 0.8 * delta_T;  // 0.8 damping factor

    // Recompute residual
    K = partition_coefficients(pressure, T_solidus, A, B, L, T0, R);
    residual = compute_residual(composition, K, compute_solidus);

    ++n;
}

if (n == max_iterations) {
    std::cerr << "!!! Newton solver did not converge after "
              << n << " iterations. Final residual = " << residual << " !!!" << std::endl;
}

      //auto tend = std::chrono::high_resolution_clock::now();
      //std::chrono::duration<double> elapsed = tend - start;
      //totaln = totaln + n;
      //totaltime = totaltime + elapsed.count();
      //std::cout<<n<<" "<<T_solidus<<" "<<elapsed.count()<<" "<<totaln<<" "<<totaltime<<std::endl;
      return T_solidus;
    }

    template <int dim>
    std::vector<double>
    VolatilesMelt<dim>::melting_temperatures(const double pressure,
                                             std::vector<double> A,
                                             std::vector<double> B,
                                             std::vector<double> L,
                                             std::vector<double> T0) const
    {
        std::vector<double> Tm (n_components);

        const double Pmax = 6e9;
        if(use_simons_law)
        {
          for (unsigned int i=0; i<n_components; ++i)
            Tm[i]  =  T0[i] * std::pow(1.0 + pressure/A[i], 1.0/B[i]);
        }
        else
        {
          if (pressure <= Pmax)
            for (unsigned int i=0; i<n_components; ++i) 
              Tm[i]  =  T0[i] + A[i] * pressure + B[i] * pressure * pressure;
          else
          {
            // safeguard: continue melting point with linear slope above Pmax
            const double dP = 1e7;
            for (unsigned int i=0; i<n_components; ++i) 
            {
              const double T0_at_Pmax = T0[i] + A[i] * Pmax + B[i] * Pmax * Pmax;
              const double dTdP = ((A[i]*Pmax + B[i]*Pmax*Pmax) - (A[i]*(Pmax-1e7) + B[i]*(Pmax-1e7)*(Pmax-1e7)))/dP;
              Tm[i] = T0_at_Pmax + dTdP * (pressure-Pmax);
            }
          }
        }

        return Tm;
    }


    template <int dim>
    std::vector<double>
    VolatilesMelt<dim>::partition_coefficients(const double pressure, 
                                                const double temperature,
                                                std::vector<double> A,
                                                std::vector<double> B,
                                                std::vector<double> L,
                                                std::vector<double> T0,
                                                std::vector<double> R) const
    {
        std::vector<double> K (n_components);

        std::vector<double> Tm = melting_temperatures(pressure, A, B, L, T0);

        // Parameterization after Rudge, Bercovici, & Spiegelman (2010)
        for (unsigned int i=0; i<n_components; ++i) 
        {
            // L/T gives the change in entropy.
            double Ls = L[i]/T0[i]*temperature;
            K[i] = std::exp(Ls/R[i] * (1./temperature - 1./Tm[i]));
        }

        return K;
    }

    template <int dim>
    double 
    VolatilesMelt<dim>::compute_residual (std::vector<double> composition,
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
    double
    VolatilesMelt<dim>::
    get_reaction_rate (const double old_value,
                                     const double reaction_rate,
                                     const double time_step,
                                     const double depth,
                                     const double x) const
    {
      // Ensure composition stays between 1 and 0.
      double final_reaction_rate = reaction_rate;

      if (old_value > 1.0 || old_value + reaction_rate * time_step > 1.0)
        final_reaction_rate = (1.0 - old_value) / time_step;
      else if (old_value < 0.0 || old_value + reaction_rate * time_step < 0.0)
        final_reaction_rate = -old_value / time_step;

      if(use_extraction_patch)
        if(depth <  extraction_depth && x < extraction_width)
              final_reaction_rate = 0.0;

      return final_reaction_rate;
    }

    template <int dim>
    double
    VolatilesMelt<dim>::linear_interpolation(double P,
                                                 double P1,
                                                 double P2,
                                                 double v1,
                                                 double v2) const
    {
      // Linear interpolation
      double t = (P - P1) / (P2 - P1);
      double value = v1 + t * (v2 - v1);

      // Clamp so that we don't extrapolate beyond val1–val2
      return std::clamp(value, std::min(v1, v2), std::max(v1, v2));

    }

      template <int dim>
      void
      VolatilesMelt<dim>::initialize_component_values (std::vector<double> &Ac,
                                                          std::vector<double> &Bc,
                                                          std::vector<double> &Lc,
                                                          std::vector<double> &Tc,
                                                          std::vector<double> &Rc,
                                                          const double pressure) const
      {
        double P1 = 1.8e9; double P2 = 2.1e9;

        // High pressure parameters for cmorb
        double T2 = 1238; double A2 = 10e-9; double B2=4e-18;
        for (unsigned int i=0; i<n_components; ++i) 
        {
            Ac[i] = Ai[i];
            Bc[i] = Bi[i];
            Lc[i] = Li[i];
            Tc[i] = T0i[i];
            Rc[i] = Ri[i];
        }

      /*if(pressure >= P2)
       {
          Ac[2] = A2;
          Bc[2] = B2;
          Tc[2] = T2;
       }
       else if (pressure < P2 && pressure > P1)
       {
          Ac[2] = linear_interpolation(pressure, P1, P2, Ai[2], A2);
          Bc[2] = linear_interpolation(pressure, P1, P2, Bi[2], B2);
          Tc[2] = linear_interpolation(pressure, P1, P2, T0i[2], T2);
       }*/

      }        


      template <int dim>
      void
      VolatilesMelt<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Reaction model");
          {
          prm.enter_subsection("Volatiles model");
          {
              prm.declare_entry ("T0", "1085.7",
                                 Patterns::List(Patterns::Double (0.)),
                                 "Constant parameter in the quadratic "
                                 "function that approximates the solidus "
                                 "of peridotite. "
                                 "Units: \\si{\\degreeCelsius}.");
              prm.declare_entry ("A", "1.329e-7",
                                 Patterns::List(Patterns::Double ()),
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

              prm.declare_entry ("Fluid density difference", "500",
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
              prm.declare_entry ("Maximum permeability", "1e-1",
                                Patterns::Double(),
                                "Maximum permeability of the solid host rock."
                                "Units: \\si{\\meter\\squared}.");
              prm.declare_entry ("Maximum pressure for melt", "4.75e12",
                                Patterns::Double (0.),
                                "The value of the constant melt viscosity $\\viscosity_fluid$. Units: \\si{\\pascal\\second}.");
              prm.declare_entry ("Extraction depth", "4000",
                                Patterns::Double (),
                                "The value of the constant melt viscosity $\\viscosity_fluid$. Units: \\si{\\pascal\\second}.");
              prm.declare_entry ("Extraction width", "4000",
                                Patterns::Double (),
                                "The value of the constant melt viscosity $\\viscosity_fluid$. Units: \\si{\\pascal\\second}.");
              prm.declare_entry ("Compaction to shear viscosity ratio", "10",
                                Patterns::Double (),
                                "The value of the constant melt viscosity $\\viscosity_fluid$. Units: \\si{\\pascal\\second}.");
              prm.declare_entry ("Use fractional melting", "false",
                             Patterns::Bool (),
                             "Whether to use fractional or batch melting.");
              prm.declare_entry ("Use extraction patch", "false",
                             Patterns::Bool (),
                             "Whether to use fractional or batch melting.");
              prm.declare_entry ("Use simons law", "false",
                             Patterns::Bool (),
                             "Whether to use fractional or batch melting.");
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
      
      


      template <int dim>
      void
      VolatilesMelt<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Reaction model");
        {
        prm.enter_subsection("Volatiles model");
        {
            n_components = prm.get_integer("Number of components");
            T0i = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("T0"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            Ai = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("A"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");                                                              
            Bi = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("B"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            Li = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("L"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            Ri = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("R"))),
                                                                          n_components,
                                                                          "Thermal diffusivities");
            fluid_density_difference         = prm.get_double ("Fluid density difference");
            extraction_depth        = prm.get_double ("Extraction depth");
            extraction_width        = prm.get_double ("Extraction width");
            pressure_max         = prm.get_double ("Maximum pressure for melt");
            compaction_viscosity_ratio         = prm.get_double ("Compaction to shear viscosity ratio");
            melting_time_scale         = prm.get_double ("Melting time scale for operator splitting");
            use_fractional_melting        = prm.get_bool ("Use fractional melting");
            use_extraction_patch        = prm.get_bool ("Use extraction patch");
            use_simons_law     = prm.get_bool ("Use simons law");

            if (this->convert_output_to_years() == true)
              melting_time_scale *= year_in_seconds;

            // Get ID for all fields. Should I do this once here or do it in the main function?
            AssertThrow(this->introspection().compositional_name_exists("porosity"), ExcMessage("A porosity field is needed to use the co2 plugin."));
            porosity_idx = this->introspection().compositional_index_for_name("porosity");

            AssertThrow(this->introspection().compositional_name_exists("morb_cl"), ExcMessage("A morb_cl field is needed to use the co2 plugin."));
            mcl_idx = this->introspection().compositional_index_for_name("morb_cl");

            AssertThrow(this->introspection().compositional_name_exists("morb_cs"), ExcMessage("A morb_cs field is needed to use the co2 plugin."));
            mcs_idx = this->introspection().compositional_index_for_name("morb_cs");
          
            AssertThrow(this->introspection().compositional_name_exists("cmorb_cl"), ExcMessage("A cmorb_cl field is needed to use the co2 plugin."));
            ccl_idx = this->introspection().compositional_index_for_name("cmorb_cl");

            AssertThrow(this->introspection().compositional_name_exists("cmorb_cs"), ExcMessage("A cmorb_cs field is needed to use the co2 plugin."));
            ccs_idx = this->introspection().compositional_index_for_name("cmorb_cs");

            AssertThrow(this->introspection().compositional_name_exists("hmorb_cl"), ExcMessage("A hmorb_cl field is needed to use the co2 plugin."));
            hcl_idx = this->introspection().compositional_index_for_name("hmorb_cl");

            AssertThrow(this->introspection().compositional_name_exists("hmorb_cs"), ExcMessage("A hmorb_cs field is needed to use the co2 plugin."));
            hcs_idx = this->introspection().compositional_index_for_name("hmorb_cs");

            xi_0                       = prm.get_double ("Reference bulk viscosity");
            viscosity_fluid            = prm.get_double ("Reference melt viscosity");
            thermal_bulk_viscosity_exponent = prm.get_double ("Thermal bulk viscosity exponent");
            alpha_phi                  = prm.get_double ("Exponential melt weakening factor");
            melt_compressibility       = prm.get_double ("Melt compressibility");
            melt_bulk_modulus_derivative = prm.get_double ("Melt bulk modulus derivative");
            reference_permeability     = prm.get_double ("Reference permeability");
            maximum_permeability     = prm.get_double ("Maximum permeability");
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
    template class VolatilesMelt<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
