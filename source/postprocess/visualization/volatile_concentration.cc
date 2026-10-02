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

#include <aspect/postprocess/visualization/volatile_concentration.h>
#include <aspect/melt.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/simulator.h>
#include <aspect/utilities.h>
#include <aspect/material_model/interface.h>

#include <deal.II/numerics/data_out.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      VolatilesConcentration<dim>::
      VolatilesConcentration ()
        :
        DataPostprocessor<dim> (),
        Interface<dim>("")  // a fraction, so physical units "1"
      {}

      template <int dim>
      std::vector<std::string>
      VolatilesConcentration<dim>::
      get_names () const
      {
        std::vector<std::string> solution_names;
        solution_names.emplace_back("co2_ppm");
        solution_names.emplace_back("co2_mass");
        solution_names.emplace_back("h20_ppm");
        solution_names.emplace_back("h20_mass");
        solution_names.emplace_back("morb_ppm");
        solution_names.emplace_back("morb_mass");
        solution_names.emplace_back("div_u_phys");
        solution_names.emplace_back("div_u_solution");
        return solution_names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      VolatilesConcentration<dim>::
      get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation(8,
            DataComponentInterpretation::component_is_scalar);

        return interpretation;
      }


      template <int dim>
      UpdateFlags
      VolatilesConcentration<dim>::
      get_needed_update_flags () const
      {
        return update_gradients | update_values  | update_quadrature_points;
      }

      template <int dim>
      void
      VolatilesConcentration<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double>> &computed_quantities) const
      {
        AssertThrow(this->include_melt_transport()==true,
                    ExcMessage("'Include melt transport' has to be on when using melt transport postprocessors."));

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());

        // Set use_strain_rates to true since the compaction viscosity might also depend on the strain rate.
        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points, this->n_compositional_fields());
        MeltHandler<dim>::create_material_model_outputs(out);

        this->get_material_model().evaluate(in, out);
        const std::shared_ptr<MaterialModel::MeltOutputs<dim>> fluid_out
          = out.template get_additional_output_object<MaterialModel::MeltOutputs<dim>>();
        AssertThrow(fluid_out != nullptr,
                    ExcMessage("Need MeltOutputs from the material model for computing the melt properties."));

        const double p_c_scale = Plugins::get_plugin_as_type<const MaterialModel::MeltInterface<dim>>(this->get_material_model()).p_c_scale(in,
                                 out,
                                 this->get_melt_handler(),
                                 true);

        // Fetch compaction pressure at all q-points (component from introspection).
         const unsigned int p_c_component = this->introspection().variable("compaction pressure").first_component_index;
        // Cell-averaged, melt-cell-gated divergence, matching melt.cc.
        double divergence_u = 0.0;
        double div_u2 = 0.0;
        const bool is_melt_cell = (p_c_scale > 0.0);

        if (is_melt_cell)
          for (unsigned int q=0; q<n_quadrature_points; ++q)
            {
              const double xi = fluid_out->compaction_viscosities[q];
              const double p_c = input_data.solution_values[q][p_c_component];
              const double div_u_phys = (xi > 0.0 ? - p_c_scale * p_c / xi : 0.0);
              divergence_u += div_u_phys * 1./n_quadrature_points;

            double div_at_q = 0.0;
            for (unsigned int d=0; d<dim; ++d)
              {
                const unsigned int u_component = this->introspection().component_indices.velocities[d];
                div_at_q += input_data.solution_gradients[q][u_component][d];
              }
            div_u2 += div_at_q * 1./n_quadrature_points;
            }

          for (unsigned int q=0; q<n_quadrature_points; ++q)
            {
              std::vector<double> composition(this->n_compositional_fields());
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                composition[c] = input_data.solution_values[q][this->introspection().component_indices.compositional_fields[c]];

              // Get the compositional values, and limit between 0 and 1.
              const double cmorb_cl =  std::max(0.0, std::min(composition[ccl_idx],1.0));
              const double cmorb_cs =  std::max(0.0, std::min(composition[ccs_idx],1.0));
              const double hmorb_cl =  std::max(0.0, std::min(composition[hcl_idx],1.0));
              const double hmorb_cs =  std::max(0.0, std::min(composition[hcs_idx],1.0));
              const double morb_cl =  std::max(0.0, std::min(composition[mcl_idx],1.0));
              const double morb_cs =  std::max(0.0, std::min(composition[mcs_idx],1.0));
              const double Fvol =  std::max(0.0, std::min(composition[porosity_idx],1.0));
              double rho_s = out.densities[q];
              const double rho_l = fluid_out->fluid_densities[q];

              // We track the volume fraction of melt, convert to mass fraction here.
              const double avg_rho = Fvol*rho_l + (1 - Fvol)*rho_s;
              const double Fmass = Fvol*rho_l/avg_rho;     
              
              // request update_gradients in get_needed_update_flags(), then:
              //

              // Compute ppm of different compositions. Here we use the C_bar calculated from the mass fraction,
              // and multiply it by the weight percent that is co2 or h2o, and then apply a scaling factor.
              // Generally, compositions would be given in a percentage, as well as e.g., cwt, these would both need to
              // be divided by 100 and then everything multiplied by 1e6 to scale to ppm. However, in our case 
              // compositions are already not given in the percent value, so we can remove that from the scaling factor.
              computed_quantities[q](0) = (Fmass * cmorb_cl + (1 - Fmass)*cmorb_cs) * cwt/100 * 1e6;

              // Next, we output the mass of the composition. In this case, if we want to use the real
              // densities we multiply by the volume fraction of melt. If we wanted to use the mass fraction,
              // we would multiply the whole thing by the avg_rho instead.
              computed_quantities[q](1) = (Fvol * cmorb_cl * rho_l + (1 - Fvol)*cmorb_cs*rho_s) * cwt/100;

              // H2O content.
              computed_quantities[q](2) = (Fmass * hmorb_cl + (1 - Fmass) * hmorb_cs) * hwt/100 * 1e6;
              computed_quantities[q](3) = (Fvol * hmorb_cl * rho_l + (1 - Fvol)*hmorb_cs*rho_s) * hwt/100;

              // MORB content.
              computed_quantities[q](4) = (Fmass * morb_cl + (1 - Fmass)*morb_cs) * 100/100 * 1e6;
              computed_quantities[q](5) = (Fvol * morb_cl * rho_l + (1 - Fvol)*morb_cs*rho_s) * 100/100;


              const double xi = fluid_out->compaction_viscosities[q];
              const double p_c = input_data.solution_values[q][p_c_component];
              const double div_u_phys = (xi > 0.0 ? - p_c_scale * p_c / xi : 0.0);
              computed_quantities[q](6) = div_u_phys;
              computed_quantities[q](7) = div_u2;


        }
      }

   

      template <int dim>
      void
      VolatilesConcentration<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Volatiles concentration");
            {
              prm.declare_entry ("Co2 weight percent", "20",
                                 Patterns::List(Patterns::Double (0.)),
                                 "Constant parameter in the quadratic "
                                 "function that approximates the solidus "
                                 "of peridotite. "
                                 "Units: \\si{\\degreeCelsius}.");
              prm.declare_entry ("H2o weight percent", "5",
                                 Patterns::List(Patterns::Double (0.)),
                                 "Prefactor of the linear pressure term "
                                 "in the quadratic function that approximates "
                                 "the solidus of peridotite. "
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
              prm.declare_entry ("Number of components", "2",
                                 Patterns::Integer(),
                                 "Prefactor of the linear pressure term "
                                 "in the quadratic function that approximates "
                                 "the  lherzolite liquidus used for "
                                 "calculating the fraction of peridotite-"
                                 "derived melt. "
                                 "\\si{\\degreeCelsius\\per\\pascal}.");  

            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      VolatilesConcentration<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Visualization");
          {
            prm.enter_subsection("Volatiles concentration");
            {
            fluid_density_difference         = prm.get_double ("Fluid density difference");
            rho_s         = prm.get_double ("Solid density");
            cwt         = prm.get_double ("Co2 weight percent");
            hwt         = prm.get_double ("H2o weight percent");
            n_components = prm.get_integer("Number of components");

            // Get ID for all fields. Should I do this once here or do it in the main function?
            AssertThrow(this->introspection().compositional_name_exists("porosity"), ExcMessage("A porosity field is needed to use the co2 plugin."));
            porosity_idx = this->introspection().compositional_index_for_name("porosity");
          
            AssertThrow(this->introspection().compositional_name_exists("cmorb_cl"), ExcMessage("A cmorb_cl field is needed to use the co2 plugin."));
            ccl_idx = this->introspection().compositional_index_for_name("cmorb_cl");

            AssertThrow(this->introspection().compositional_name_exists("cmorb_cs"), ExcMessage("A cmorb_cs field is needed to use the co2 plugin."));
            ccs_idx = this->introspection().compositional_index_for_name("cmorb_cs");

            AssertThrow(this->introspection().compositional_name_exists("hmorb_cl"), ExcMessage("A hmorb_cl field is needed to use the co2 plugin."));
            hcl_idx = this->introspection().compositional_index_for_name("hmorb_cl");

            AssertThrow(this->introspection().compositional_name_exists("hmorb_cs"), ExcMessage("A hmorb_cs field is needed to use the co2 plugin."));
            hcs_idx = this->introspection().compositional_index_for_name("hmorb_cs");

            AssertThrow(this->introspection().compositional_name_exists("morb_cl"), ExcMessage("A hmorb_cl field is needed to use the co2 plugin."));
            mcl_idx = this->introspection().compositional_index_for_name("morb_cl");

            AssertThrow(this->introspection().compositional_name_exists("morb_cs"), ExcMessage("A hmorb_cs field is needed to use the co2 plugin."));
            mcs_idx = this->introspection().compositional_index_for_name("morb_cs");

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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(VolatilesConcentration,
                                                  "volatile concentrations", // TODO write down equations here
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
