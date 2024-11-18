/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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


#include </home/bbpdneu1/software/aspect/aspect/melt_post_plugin/melt_components.h>
#include <aspect/melt.h>

#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      MeltComponents<dim>::
      MeltComponents ()
        :
        DataPostprocessor<dim> (),
        Interface<dim>("")  // a fraction, so physical units "1"
      {}

      template <int dim>
      std::vector<std::string>
      MeltComponents<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;
        solution_names.emplace_back("Melt");
        solution_names.emplace_back("mCl");
        solution_names.emplace_back("mCs");
        solution_names.emplace_back("cCl");
        solution_names.emplace_back("cCs");
        return solution_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      MeltComponents<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(5,
            DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }


      template <int dim>
      UpdateFlags
      MeltComponents<dim>::
      get_needed_update_flags () const
      {
        return update_values | update_quadrature_points;
      }

      template <int dim>
      void
      MeltComponents<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        //Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

          for (unsigned int q=0; q<n_quadrature_points; ++q)
            {
              const double pressure    = input_data.solution_values[q][this->introspection().component_indices.pressure];
              const double temperature = input_data.solution_values[q][this->introspection().component_indices.temperature];
              std::vector<double> composition(this->n_compositional_fields());
              std::vector<double> avg_rho(n_components);
              std::vector<double> F(n_components);
              std::vector<double> C_bar(n_components);

              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                composition[c] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[c]];

              //std::max(0.0, std::min(1.0, C_bar[1] / (feq + (1 - feq) * K[1])));
              // Get the compositional values, and set to zero if they don't exist.
              double morb_cl =  this->introspection().compositional_name_exists("morb_cl")
                                ?
                                std::max(0.0, std::min(composition[this->introspection().compositional_index_for_name("morb_cl")],1.0))
                                :
                                (n_components==3 ? 0. : -1);
              double morb_cs =  this->introspection().compositional_name_exists("morb_cs")
                                ?
                                std::max(0.0, std::min(composition[this->introspection().compositional_index_for_name("morb_cs")],1.0))
                                :
                                (n_components==3 ? 0. : -1);
              double cmorb_cl = this->introspection().compositional_name_exists("cmorb_cl")
                                ?
                                std::max(0.0, std::min(composition[this->introspection().compositional_index_for_name("cmorb_cl")],1.0))
                                :
                                (n_components==3 ? 0. : -1);
              double cmorb_cs = this->introspection().compositional_name_exists("cmorb_cs")
                                ?
                                std::max(0.0, std::min(composition[this->introspection().compositional_index_for_name("cmorb_cs")],1.0))
                                :
                                (n_components==3 ? 0. : -1); 
              double F_int = this->introspection().compositional_name_exists("feq")
                                ?
                                std::max(0.0, std::min(composition[this->introspection().compositional_index_for_name("feq")],1.0))
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
              std::vector<double> order_comp;
              for (unsigned int i=0; i<n_components; ++i)
              {
                // May use this if convert to volume fraction, but paper uses mass.
                // Note: Paper does not have different densities for component,
                // There is one liquid density and one solid density.
                //avg_rho[i] = porosity*rho_l[i] + (1 - porosity)*rho_s[i];
                C_bar[i] = F_int*c_l[i] + (1-F_int)*c_s[i];
                order_comp.push_back(C_bar[i]);
              } 
           
              const double T_solidus = T_solidus_liquidus(pressure, order_comp, true);
              const double T_liquidus = T_solidus_liquidus(pressure, order_comp, false);

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

              if(temperature <= T_solidus)
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
              feq = std::max(0.0, std::min(1.0, feq));
              double cl = std::max(0.0, std::min(1.0, C_bar[1] / (feq + (1 - feq) * K[1])));
              double cs = std::max(0.0, std::min(1.0, C_bar[1] / (feq / K[1] + (1 - feq))));
              double cl2 = std::max(0.0, std::min(1.0, C_bar[2] / (feq + (1 - feq) * K[2])));
              double cs2 = std::max(0.0, std::min(1.0, C_bar[2] / (feq / K[2] + (1 - feq))));

              const unsigned int melt_idx = this->introspection().compositional_index_for_name("feq");
              const unsigned int cl_idx = this->introspection().compositional_index_for_name("morb_cl");
              const unsigned int cs_idx = this->introspection().compositional_index_for_name("morb_cs");
              const unsigned int cl2_idx = this->introspection().compositional_index_for_name("cmorb_cl");
              const unsigned int cs2_idx = this->introspection().compositional_index_for_name("cmorb_cs");

              //std::cout<<C_bar[1]<<" "<<feq<<" "<<K[1]<<" "<<cs<<" "<<cl<<std::endl;
              computed_quantities[q](0) = feq;
              computed_quantities[q](1) = cl;
              computed_quantities[q](2) = cs;
              computed_quantities[q](3) = cl2;
              computed_quantities[q](4) = cs2;
              
              // WHAT TO OUTPUT: convert f back to volume fraction

            }
        //double T = Tm;

        //double r = ((C*K)^2)^-1;

        //return 0.0;

      }

    template <int dim>
    std::vector<double>
    MeltComponents<dim>::melting_temperatures(const double pressure) const
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
    MeltComponents<dim>::partition_coefficients(const double pressure, 
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
    MeltComponents<dim>::compute_residual (std::vector<double> composition,
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
      MeltComponents<dim>::
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
        const unsigned int max_iterations =  500;   // maximum number of iterations
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
            std::cerr << "!!! Newton solver for solidus/liquidus T has not converged after " << max_iterations << " iterations !!!" << std::endl;
            break;
          }
        }
      return T_solidus;
    } 

      template <int dim>
      void
      MeltComponents<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt components");
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

              prm.declare_entry ("Porosity", "0",
                    Patterns::Double(0,1),
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
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      MeltComponents<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Melt components");
            {

            n_components = prm.get_integer("Number of components");
            porosity = prm.get_double("Porosity");
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
            rho_l = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Fluid density"))),
                                                              n_components,
                                                              "Thermal diffusivities");

            rho_s = Utilities::possibly_extend_from_1_to_N (Utilities::string_to_double(Utilities::split_string_list(prm.get("Solid density"))),
                                                              n_components,
                                                              "Thermal diffusivities");
            }
            prm.leave_subsection();
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
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(MeltComponents,
                                                  "melt components", // TODO write down equations here
                                                  "A visualization output object that generates output "
                                                  "for the melt fraction at the temperature and "
                                                  "pressure of the current point. If the material model computes "
                                                  "a melt fraction, this is the quantity that will be visualized. "
                                                  "Otherwise, a specific parametrization for batch melting "
                                                  "(as described in the following) will be used. "
                                                  "It does not take into account latent heat. "
                                                  "If there are no compositional fields, or no fields called 'pyroxenite', "
                                                  " this postprocessor will visualize the melt fraction of peridotite "
                                                  "(calculated using the anhydrous model of Katz, 2003). "
                                                  "If there is a compositional field called 'pyroxenite', the "
                                                  "postprocessor assumes that this compositional "
                                                  "field is the content of pyroxenite, and will visualize "
                                                  "the melt fraction for a mixture of peridotite and pyroxenite "
                                                  "(using the melting model of Sobolev, 2011 for pyroxenite). "
                                                  "All the parameters that were used in these calculations "
                                                  "can be changed in the input file, the most relevant maybe "
                                                  "being the mass fraction of Cpx in peridotite in the Katz "
                                                  "melting model (Mass fraction cpx), which right now has a "
                                                  "default of 15\\%. "
                                                  "The corresponding $p$-$T$-diagrams can be generated by running "
                                                  "the tests melt\\_postprocessor\\_peridotite and "
                                                  "melt\\_postprocessor\\_pyroxenite."
                                                  "\n\n"
                                                  "Physical units: None.")
    }
  }
}
