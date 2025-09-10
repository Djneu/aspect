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

#include </mnt/vast-nhr/home/derekjohn.neuharth/u16318/software/co2/aspect/melt_volatiles_plugins/volatiles_melt_postprocess.h>
#include <aspect/melt.h>
#include <deal.II/base/parameter_handler.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VolatilesMeltPost<dim>::
      VolatilesMeltPost ()
        :
        DataPostprocessor<dim> (),
        Interface<dim>("")  // a fraction, so physical units "1"
      {}

      template <int dim>
      std::vector<std::string>
      VolatilesMeltPost<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;
        solution_names.emplace_back("Melt");
        solution_names.emplace_back("mCl");
        solution_names.emplace_back("mCs");
        solution_names.emplace_back("cCl");
        solution_names.emplace_back("cCs");
        solution_names.emplace_back("hCl");
        solution_names.emplace_back("hCs");
        solution_names.emplace_back("cppm");
        solution_names.emplace_back("hppm");
        return solution_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      VolatilesMeltPost<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(9,
            DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }


      template <int dim>
      UpdateFlags
      VolatilesMeltPost<dim>::
      get_needed_update_flags () const
      {
        return update_values | update_quadrature_points;
      }

      template <int dim>
      void
      VolatilesMeltPost<dim>::
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
              std::vector<double> C_bar(n_components);

              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                composition[c] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[c]];


              // Get the compositional values, and limit between 0 and 1.
              double morb_cl =  std::max(0.0, std::min(composition[mcl_idx],1.0));
              double morb_cs =  std::max(0.0, std::min(composition[mcs_idx],1.0));
              double cmorb_cl =  std::max(0.0, std::min(composition[ccl_idx],1.0));
              double cmorb_cs =  std::max(0.0, std::min(composition[ccs_idx],1.0));
              double hmorb_cl =  std::max(0.0, std::min(composition[hcl_idx],1.0));
              double hmorb_cs =  std::max(0.0, std::min(composition[hcs_idx],1.0));
              double Fvol_old =  0.0; //std::max(0.0, std::min(composition[melt_idx],1.0));
              const double rho_l = rho_s - fluid_density_difference;                                                                           

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
   
        // Now that things are ordered, find the bulk composition for each component.
        for (unsigned int i=0; i<n_components; ++i)
            C_bar[i] = Fmass_old*c_l[i] + (1-Fmass_old)*c_s[i];
            
        //std::cout<<C_bar[0]<<" "<<C_bar[1]<<" "<<C_bar[2]" "<<C_bar[0]<<std::endl;
        // Define parameters that will be returned.
        double Fmass_new = 0.0;
           
      // Calculate the equilibrium values and reaction rates
      // if we are below the maximum solidus pressure.
         double mcl = 0.0;
        double ccl = 0.0;    
        double hcl = 0.0;
        double mcs = 0.0;
        double ccs = 0.0;
        double hcs = 0.0;
        const double T_solidus = T_solidus_liquidus(pressure, C_bar, true);
        const double T_liquidus = T_solidus_liquidus(pressure, C_bar, false);

          std::vector<double> Tm = melting_temperatures(pressure);
          std::vector<double> K = partition_coefficients(pressure, std::max(T_solidus,std::min(T_liquidus,temperature)));
          
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

          // Liquid values
          double dcl = std::max(0.0, std::min(1.0, C_bar[0] / (Fmass_new + (1 - Fmass_new) * K[0])));
          mcl = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new + (1 - Fmass_new) * K[1])));
          ccl = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new + (1 - Fmass_new) * K[2])));     
          hcl = std::max(0.0, std::min(1.0, C_bar[3] / (Fmass_new + (1 - Fmass_new) * K[3])));
          
          // Solid values
          double dcs = std::max(0.0, std::min(1.0, C_bar[0] / (Fmass_new / K[0] + (1 - Fmass_new))));
          mcs = std::max(0.0, std::min(1.0, C_bar[1] / (Fmass_new / K[1] + (1 - Fmass_new))));
          ccs = std::max(0.0, std::min(1.0, C_bar[2] / (Fmass_new / K[2] + (1 - Fmass_new))));
          hcs = std::max(0.0, std::min(1.0, C_bar[3] / (Fmass_new / K[3] + (1 - Fmass_new))));

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

          computed_quantities[q](0) = Fmass_new; //Fmass_new*(rho_l/avg_rho);
          computed_quantities[q](1) = mcl;
          computed_quantities[q](2) = mcs;
          computed_quantities[q](3) = ccl;
          computed_quantities[q](4) = ccs;
          computed_quantities[q](5) = hcl;
          computed_quantities[q](6) = hcs;   
          computed_quantities[q](7) = (Fmass_new * ccl + (1 - Fmass_new)*ccs) * 20/100 * 1e6;
          computed_quantities[q](8) = (Fmass_new * hcl + (1 - Fmass_new)*hcs) * 5/100 * 1e6;;
        }
      }

    template <int dim>
    std::vector<double>
    VolatilesMeltPost<dim>::melting_temperatures(const double pressure) const
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
    VolatilesMeltPost<dim>::partition_coefficients(const double pressure, 
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
        {
            double Ls = L[i]/T0[i]*temperature;
            K[i] = std::exp(Ls/R[i] * (1./(temperature) - 1./(Tm[i])));
        }

        return K;
    }

    template <int dim>
    double 
    VolatilesMeltPost<dim>::compute_residual (std::vector<double> composition,
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
      VolatilesMeltPost<dim>::
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
    K = partition_coefficients(pressure, T_solidus + eps_T);
    double residual_plus = compute_residual(composition, K, compute_solidus);

    K = partition_coefficients(pressure, T_solidus - eps_T);
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
    K = partition_coefficients(pressure, T_solidus);
    residual = compute_residual(composition, K, compute_solidus);

    ++n;
}

if (n == max_iterations) {
    std::cerr << "??? Newton solver did not converge after "
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
      void
      VolatilesMeltPost<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Volatiles melt");
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

              prm.declare_entry ("Fluid density difference", "500",
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

              prm.declare_entry ("Maximum pressure for melt", "4.75e9",
                                Patterns::Double (0.),
                                "The value of the constant melt viscosity $\\viscosity_fluid$. Units: \\si{\\pascal\\second}.");

              prm.declare_entry ("Use simons law", "false",
                             Patterns::Bool (),
                             "Whether to use fractional or batch melting.");

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      VolatilesMeltPost<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Volatiles melt");
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
            fluid_density_difference         = prm.get_double ("Fluid density difference");
            rho_s         = prm.get_double ("Solid density");
            pressure_max         = prm.get_double ("Maximum pressure for melt");
            use_simons_law     = prm.get_bool ("Use simons law");

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

            AssertThrow(this->introspection().compositional_name_exists("hmorb_cl"), ExcMessage("A hmorb_cl field is needed to use the co2 plugin."));
            hcl_idx = this->introspection().compositional_index_for_name("hmorb_cl");

            AssertThrow(this->introspection().compositional_name_exists("hmorb_cs"), ExcMessage("A hmorb_cs field is needed to use the co2 plugin."));
            hcs_idx = this->introspection().compositional_index_for_name("hmorb_cs");

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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VolatilesMeltPost,
                                                  "volatiles melt post", // TODO write down equations here
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
