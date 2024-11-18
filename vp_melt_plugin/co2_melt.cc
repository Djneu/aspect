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
      void
      Co2Melt<dim>::
      calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                      typename Interface<dim>::MaterialModelOutputs &out) const
      {
          ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();

          for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
            {

              // Define component dependent parameters
              std::vector<double> F(n_components);  
              std::vector<double> C_bar(n_components);

              std::vector<double> composition(this->n_compositional_fields());
              const double temperature = in.temperature[q];
              double pressure    = in.pressure[q];

              // Set pressure to surface pressure if it is below zero.
              if(pressure < 0)
                pressure = 101325;

              // Get the compositional values, and limit between 0 and 1.
              // If they don't exist, assume it is a 2 component system and
              // set to -1 so that they are not added to component dependent parameters.
              double morb_cl =  this->introspection().compositional_name_exists("morb_cl")
                                ?
                                std::max(0.0, std::min(in.composition[q][this->introspection().compositional_index_for_name("morb_cl")],1.0))
                                :
                                (n_components==3 ? 0. : -1);
              double morb_cs =  this->introspection().compositional_name_exists("morb_cs")
                                ?
                                std::max(0.0, std::min(in.composition[q][this->introspection().compositional_index_for_name("morb_cs")],1.0))
                                :
                                (n_components==3 ? 0. : -1);
              double cmorb_cl = this->introspection().compositional_name_exists("cmorb_cl")
                                ?
                                std::max(0.0, std::min(in.composition[q][this->introspection().compositional_index_for_name("cmorb_cl")],1.0))
                                :
                                (n_components==3 ? 0. : -1);
              double cmorb_cs = this->introspection().compositional_name_exists("cmorb_cs")
                                ?
                                std::max(0.0, std::min(in.composition[q][this->introspection().compositional_index_for_name("cmorb_cs")],1.0))
                                :
                                (n_components==3 ? 0. : -1); 
              double F_int = this->introspection().compositional_name_exists("feq")
                                ?
                                std::max(0.0, std::min(in.composition[q][this->introspection().compositional_index_for_name("feq")],1.0))
                                :
                                0.;                                                                         

              // Combine liquid and solid comps, and then calculate l and s dunite.
              std::vector<double> cl_comp = {morb_cl, cmorb_cl};
              std::vector<double> cs_comp = {morb_cs, cmorb_cs};
              double dunite = 1;
              double dunite_l = 1;
              for (unsigned int i=0; i<2; ++i)
              {
                if(cs_comp[i] != -1)
                  dunite = dunite - cs_comp[i];
                if(cl_comp[i] != -1)
                  dunite_l = dunite_l - cl_comp[i];
              }

             // Order fields for finding bulk composition.
             // Depending on number of components.
             std::vector<double> c_l = {dunite_l};
             std::vector<double> c_s = {dunite};
             for (unsigned int i=0; i<n_components-1; ++i)
                {
                  if(cl_comp[i] != -1)
                    c_l.push_back(cl_comp[i]);
                  if(cs_comp[i] != -1)
                    c_s.push_back(cs_comp[i]);
                }


              // 0 = dunite, 1 = morb, 2 = cmorb. In celcius, convert.
              // Now that things are ordered, find the bulk composition for each component.
              // If we track the volume of melt, need to convert back to mass.
              double avg_rho = F_int*rho_l + (1 - F_int)*rho_s;
              F_int = F_int*avg_rho/rho_l;

              std::vector<double> order_comp;
              for (unsigned int i=0; i<n_components; ++i)
              {
                // May use this if convert to volume fraction, but paper uses mass.
                // Note: Paper does not have different densities for component,
                // There is one liquid density and one solid density.
                //avg_rho[i] = porosity*rho_l[i] + (1 - porosity)*rho_s[i];
                //F[i] = porosity/(rho_l[i]/avg_rho[i]);
                C_bar[i] = F_int*c_l[i] + (1-F_int)*c_s[i];
                order_comp.push_back(C_bar[i]);
              } 
           
              const double T_solidus = T_solidus_liquidus(pressure, order_comp, true, temperature);
              const double T_liquidus = T_solidus_liquidus(pressure, order_comp, false,temperature);

              std::vector<double> Tm = melting_temperatures(pressure);
              std::vector<double> K = partition_coefficients(pressure, std::max(T_solidus,std::min(T_liquidus,temperature)));

              // newton solver for feq     
              int n      =  0;       
              const double rnorm_tol = 1e-10;
              const int its_tol   = 100;
              double feq = 0;
              double r = 1.;
              double P_max = 4.75e9;
              for (unsigned int i=0; i<n_components; ++i)
              {
                double r1 = C_bar[i] / (feq + (1-feq)*K[i]);
                double r2 = C_bar[i] / (feq/K[i] + (1-feq));
                r +=  r1 - r2;
              }
        
              // Set melt fraction to old melt fraction.
              double f = F_int;

              if(temperature <= T_solidus || pressure > P_max)
                feq = 0;
              else if(temperature >= T_liquidus)
                feq = 1;
              else
              {
                while (abs(r) > rnorm_tol) 
                {
                  double dr_df = 0;
                  double term1 = 0;
                  double term2 = 0;
                  for (unsigned int i=0; i<n_components; ++i)
                  {
                    double numerator1 = C_bar[i] * (1.0 - K[i]);
                    double denominator1 = std::pow((feq + (1.0 - feq) * K[i]), 2);
                    term1 += numerator1 / denominator1;

                    // Second term: sum(VAR.C.*(1./PAR.K - 1)./(ff./PAR.K + (1-ff)).^2, 2)
                    double numerator2 = C_bar[i] * (1.0 / K[i] - 1.0);
                    double denominator2 = std::pow((feq / K[i] + (1.0 - feq)), 2);
                    term2 += numerator2 / denominator2;
                  }

                  dr_df = -term1+term2;

                
                  //std::cout<<r<<"  "<<f<<std::endl;
                  double a = 1;
                  while (((f - a*r/dr_df) < -1e-10) || (f - a*r/dr_df > 1-1e-10))
                  {
                    a = a/2;
                    //std::cout<<f<<"  "<<a<<"  "<<"  "<<r<<"  "<<"  "<<dr_df<<"  "<<(f - a*r/dr_df)<<std::endl;
                    if (a<1e-6)
                        AssertThrow(false, ExcMessage("a too small."));
                  }

                  f = f - a*r/dr_df;
                  feq = f;

                  r = 0;
                  for (unsigned int i=0; i<n_components; ++i)
                    r += C_bar[i]/(feq + (1-feq)*K[i]) - C_bar[i]/(feq/K[i] + (1-feq));

                  n = n+1;
                  if (n==its_tol)
                    AssertThrow(false, ExcMessage("No convergence"));

                }
              }

              // Calculate VAR.Cl and VAR.Cs, and limit all between 0 and 1.
              feq = (rho_l/avg_rho)*feq; // Convert back to volume.
              feq = std::max(0.0, std::min(1.0, feq));
              double cl = std::max(0.0, std::min(1.0, C_bar[1] / (feq + (1 - feq) * K[1])));
              double cs = std::max(0.0, std::min(1.0, C_bar[1] / (feq / K[1] + (1 - feq))));
              double cl2 = std::max(0.0, std::min(1.0, C_bar[2] / (feq + (1 - feq) * K[2])));
              double cs2 = std::max(0.0, std::min(1.0, C_bar[2] / (feq / K[2] + (1 - feq))));

              if(pressure > P_max)
              {
                cl = 0.0;
                cl2 = 0.0;
                feq = 0.0;
              }

              const unsigned int melt_idx = this->introspection().compositional_index_for_name("feq");
              const unsigned int cl_idx = this->introspection().compositional_index_for_name("morb_cl");
              const unsigned int cs_idx = this->introspection().compositional_index_for_name("morb_cs");
              const unsigned int cl2_idx = this->introspection().compositional_index_for_name("cmorb_cl");
              const unsigned int cs2_idx = this->introspection().compositional_index_for_name("cmorb_cs");

                // WHAT TO OUTPUT: convert f back to volume fraction? It is mass fraction at the moment.
                // because depletion is a volume-based, and not a mass-based property that is advected,
                // additional scaling factors on the right hand side apply

                // We use the fulls compositions here so the rate includes any deviations outside of 0 and 1.
                for (unsigned int c=0; c<in.composition[q].size(); ++c)
                  {
                    out.reaction_terms[q][c] = 0.0;
                    if (reaction_rate_out != nullptr && in.requests_property(MaterialProperties::reaction_rates) && this->get_timestep_number() > 0)
                    {
                      if (c == melt_idx)
                      {
                            double rate = (feq - in.composition[q][this->introspection().compositional_index_for_name("feq")]);
                            rate = std::max(rate, -in.composition[q][this->introspection().compositional_index_for_name("feq")]);
                            reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                      }
                      else if (c == cs_idx)
                      {
                            double rate = (cs - in.composition[q][this->introspection().compositional_index_for_name("morb_cs")]);
                            rate = std::max(rate, -in.composition[q][this->introspection().compositional_index_for_name("morb_cs")]);
                            reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                      }
                      else if (c == cl_idx)
                      {
                            double rate = (cl - in.composition[q][this->introspection().compositional_index_for_name("morb_cl")]);
                            rate = std::max(rate, -in.composition[q][this->introspection().compositional_index_for_name("morb_cl")]);
                            reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                      }
                      else if (c == cs2_idx)
                      {
                            double rate = (cs2 - in.composition[q][this->introspection().compositional_index_for_name("cmorb_cs")]);
                            rate = std::max(rate, -in.composition[q][this->introspection().compositional_index_for_name("cmorb_cs")]);
                            reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                      }
                      else if (c == cl2_idx)
                      {
                            double rate = (cl2 - in.composition[q][this->introspection().compositional_index_for_name("cmorb_cl")]);
                            rate = std::max(rate, -in.composition[q][this->introspection().compositional_index_for_name("cmorb_cl")]);
                            reaction_rate_out->reaction_rates[q][c] = rate/melting_time_scale;
                      }
                      else
                        reaction_rate_out->reaction_rates[q][c] = 0.0;
                    }
                  }
            }
      }

      template <int dim>
      double
      Co2Melt<dim>::
      T_solidus_liquidus (const double pressure, 
                          std::vector<double> composition, 
                          bool compute_solidus,
                          const double temp) const
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
          // Note the step size is set to 0.5 whereas the original r_DMC implementation uses 1 
          for (unsigned int i=0; i<n_components; ++i) 
            T_solidus = T_solidus - 0.5 * residual/dresidualdT;

          // Compute partition coefficients Ki at Tsol
          K = partition_coefficients(pressure, T_solidus);

          // Get residual at T_solidus
          residual = compute_residual(composition, K, compute_solidus);

          ++n;

          if (n == max_iterations) 
          {

                          //computed_quantities[q](0) = feq;
              const unsigned int melt_idx = this->introspection().compositional_index_for_name("feq");
              const unsigned int cl_idx = this->introspection().compositional_index_for_name("morb_cl");
              const unsigned int cs_idx = this->introspection().compositional_index_for_name("morb_cs");
              const unsigned int cl2_idx = this->introspection().compositional_index_for_name("cmorb_cl");
              const unsigned int cs2_idx = this->introspection().compositional_index_for_name("cmorb_cs");
              const unsigned int ml_idx = this->introspection().compositional_index_for_name("mantle_lithosphere");

            std::cout<< compute_solidus << " " << composition[cl_idx] << " " << composition[cs_idx] << " "<<composition[cl2_idx] <<" "<<composition[cs2_idx]<<std::endl;
            std::cout<<pressure<< " "<< temp<<" " << composition[ml_idx] << " "<<T_solidus<<std::endl;
            std::cout<<dresidualdT<< " "<<residual_minus_eps_T<<" " << residual_plus_eps_T << " "<<T_solidus<<std::endl;
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

        // Implementation in r_DMC is the following instead:
        //std::vector<double> L (n_components);
        //for (unsigned int i=0; i<n_components; ++i) 
        //  L[i] = temperature * dS[i];
        // TODO: ask Tobias Keller which we should use

        std::vector<double> Tm = melting_temperatures(pressure);

        // Parameterization after Rudge, Bercovici, & Spiegelman (2010)
        for (unsigned int i=0; i<n_components; ++i) 
            K[i] = std::exp(L[i]/R[i] * (1./temperature - 1./Tm[i]));

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

              prm.declare_entry ("Fluid density", "3200",
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
              prm.declare_entry ("Melting time scale for operator splitting", "2e2",
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
