/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include "detailed_flux.h"
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/melt.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    DetailedFluxStatistics<dim>::execute (TableHandler &statistics)
    {
      const std::string unit = (this->convert_output_to_years())
                               ?
                               "kg/yr"
                               :
                               "kg/s";
      const double in_years = (this->convert_output_to_years())
                              ?
                              year_in_seconds
                              :
                              1.0;

      // Get ID for all fields. Should I do this once here or do it in the main function?
      AssertThrow(this->introspection().compositional_name_exists("porosity"), ExcMessage("A porosity field is needed to use the co2 plugin."));
      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
    
      AssertThrow(this->introspection().compositional_name_exists("cmorb_cl"), ExcMessage("A cmorb_cl field is needed to use the co2 plugin."));
      const unsigned int ccl_idx = this->introspection().compositional_index_for_name("cmorb_cl");

      AssertThrow(this->introspection().compositional_name_exists("cmorb_cs"), ExcMessage("A cmorb_cs field is needed to use the co2 plugin."));
      const unsigned int ccs_idx = this->introspection().compositional_index_for_name("cmorb_cs");

      // create a quadrature formula based on the temperature element alone.
      const Quadrature<dim-1> &quadrature_formula = this->introspection().face_quadratures.velocities;

      FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                        this->get_fe(),
                                        quadrature_formula,
                                        update_values            | update_gradients |
                                        update_normal_vectors    |
                                        update_quadrature_points | update_JxW_values);

      std::vector<Tensor<1,dim>> fluid_velocity_values(quadrature_formula.size());    
      std::vector<std::vector<double>> composition_values (this->n_compositional_fields(),std::vector<double> (quadrature_formula.size()));

      std::map<types::boundary_id, double> local_boundary_fluxes;
      std::map<types::boundary_id, double> local_boundary_fluxes_solid;
      std::map<types::boundary_id, double> local_boundary_fluxes_fluid;

      MaterialModel::MaterialModelInputs<dim> in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_face_values.n_quadrature_points, this->n_compositional_fields());

      MeltHandler<dim>::create_material_model_outputs(out);
      in.requested_properties = MaterialModel::MaterialProperties::density;
      MaterialModel::MeltOutputs<dim> *fluid_out = out.template get_additional_output<MaterialModel::MeltOutputs<dim>>();

      // for every surface face on which it makes sense to compute a
      // mass flux and that is owned by this processor,
      // integrate the normal mass flux given by the formula
      //   j =  \rho * v * n
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f))
              {
                fe_face_values.reinit (cell, f);
                // Set use_strain_rates to false since we don't need viscosity
                in.reinit(fe_face_values, cell, this->introspection(), this->get_solution());

                this->get_material_model().evaluate(in, out);

                const FEValuesExtractors::Vector ex_u_f = this->introspection().variable("fluid velocity").extractor_vector();
                fe_face_values[ex_u_f].get_function_values (this->get_solution(), fluid_velocity_values);

                double local_normal_flux = 0;
                double local_solid_flux = 0;
                double local_fluid_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                  double cmorb_cl =  std::max(0.0, std::min(in.composition[q][ccl_idx],1.0));
                  double cmorb_cs =  std::max(0.0, std::min(in.composition[q][ccs_idx],1.0));
                  double Fvol =  std::max(0.0, std::min(in.composition[q][porosity_idx],1.0));

                  double rho_s = out.densities[q];
                  double rho_l = fluid_out->fluid_densities[q];

                  double avg_rho = Fvol*rho_l + (1 - Fvol)*rho_s;
                  double Fmass = Fvol*rho_l/avg_rho; 

                  
                  local_normal_flux += (20./100 *  
                    (
                    (Fmass * cmorb_cl * rho_l * (fluid_velocity_values[q] * fe_face_values.normal_vector(q)))
                    + ((1 - Fmass) * cmorb_cs * rho_s * (in.velocity[q] * fe_face_values.normal_vector(q)))
                    )
                    * fe_face_values.JxW(q));  

                  local_solid_flux += 0.2 * (1 - Fmass) * cmorb_cs * rho_s * (in.velocity[q] * fe_face_values.normal_vector(q)) * fe_face_values.JxW(q);
                  local_fluid_flux += 0.2 * Fmass * cmorb_cl * rho_l * (fluid_velocity_values[q] * fe_face_values.normal_vector(q)) * fe_face_values.JxW(q);                   
                  }

                const types::boundary_id boundary_indicator
                  = cell->face(f)->boundary_id();
                local_boundary_fluxes[boundary_indicator] += local_normal_flux * in_years;
                local_boundary_fluxes_fluid[boundary_indicator] += local_fluid_flux * in_years;
                local_boundary_fluxes_solid[boundary_indicator] += local_solid_flux * in_years;
              }

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_boundary_fluxes;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        const std::set<types::boundary_id>
        boundary_indicators
          = this->get_geometry_model().get_used_boundary_indicators ();
        std::vector<double> local_values;
        local_values.reserve(boundary_indicators.size());
        for (const auto p : boundary_indicators)
          local_values.push_back (local_boundary_fluxes[p]);

        // then collect contributions from all processors
        std::vector<double> global_values (local_values.size());
        Utilities::MPI::sum (local_values, this->get_mpi_communicator(), global_values);

        // and now take them apart into the global map again
        unsigned int index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_fluxes[*p] = global_values[index];
      }

      // now add all of the computed mass fluxes to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      unsigned int index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_boundary_fluxes.begin();
           p != global_boundary_fluxes.end(); ++p, ++index)
        {
          const std::string name = "Outward co2 flux through boundary with indicator "
                                   + Utilities::int_to_string(p->first)
                                   + aspect::Utilities::parenthesize_if_nonempty(this->get_geometry_model()
                                                                                 .translate_id_to_symbol_name (p->first))
                                   + " (" + unit + ")";
          statistics.add_value (name, p->second);

         
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "left")
          left_total_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "right")
          right_total_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "bottom")
          bottom_total_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "top")
          top_total_flux += p->second * this->get_timestep() / year_in_seconds;


          // also make sure that the other columns filled by this object
          // all show up with sufficient accuracy and in scientific notation
          statistics.set_precision (name, 8);
          statistics.set_scientific (name, true);

          // finally have something for the screen
          screen_text.precision(4);
          screen_text << p->second << ' ' << unit
                      << (index == global_boundary_fluxes.size()-1 ? "" : ", ");
        }



      // Solid portion
      // now communicate to get the global values
      std::map<types::boundary_id, double> global_boundary_fluxes_solid;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        const std::set<types::boundary_id>
        boundary_indicators
          = this->get_geometry_model().get_used_boundary_indicators ();
        std::vector<double> local_values_solid;
        local_values_solid.reserve(boundary_indicators.size());
        for (const auto p : boundary_indicators)
          local_values_solid.push_back (local_boundary_fluxes_solid[p]);

        // then collect contributions from all processors
        std::vector<double> global_values_solid (local_values_solid.size());
        Utilities::MPI::sum (local_values_solid, this->get_mpi_communicator(), global_values_solid);

        // and now take them apart into the global map again
        index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_fluxes_solid[*p] = global_values_solid[index];
      }

      // now add all of the computed mass fluxes to the statistics object
      // and create a single string that can be output to the screen
      index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_boundary_fluxes_solid.begin();
           p != global_boundary_fluxes_solid.end(); ++p, ++index)
        {  
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "left")
          left_solid_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "right")
          right_solid_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "bottom")
          bottom_solid_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "top")
          top_solid_flux += p->second * this->get_timestep() / year_in_seconds;
        }



      // Fluid portion
      // now communicate to get the global values
      std::map<types::boundary_id, double> global_boundary_fluxes_fluid;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        const std::set<types::boundary_id>
        boundary_indicators
          = this->get_geometry_model().get_used_boundary_indicators ();
        std::vector<double> local_values_fluid;
        local_values_fluid.reserve(boundary_indicators.size());
        for (const auto p : boundary_indicators)
          local_values_fluid.push_back (local_boundary_fluxes_fluid[p]);

        // then collect contributions from all processors
        std::vector<double> global_values_fluid (local_values_fluid.size());
        Utilities::MPI::sum (local_values_fluid, this->get_mpi_communicator(), global_values_fluid);

        // and now take them apart into the global map again
        index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
          global_boundary_fluxes_fluid[*p] = global_values_fluid[index];
      }

      // now add all of the computed mass fluxes to the statistics object
      // and create a single string that can be output to the screen
      index = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_boundary_fluxes_fluid.begin();
           p != global_boundary_fluxes_fluid.end(); ++p, ++index)
        {  
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "left")
          left_fluid_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "right")
          right_fluid_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "bottom")
          bottom_fluid_flux += p->second * this->get_timestep() / year_in_seconds;
         if(this->get_geometry_model().translate_id_to_symbol_name (p->first) == "top")
          top_fluid_flux += p->second * this->get_timestep() / year_in_seconds;
        }

      statistics.add_value("Left co2 mass flow",left_total_flux);
      statistics.set_precision ("Left co2 mass flow", 7);
      statistics.set_scientific ("Left co2 mass flow", true);

      statistics.add_value("Right co2 mass flow",right_total_flux);
      statistics.set_precision ("Right co2 mass flow", 7);
      statistics.set_scientific ("Right co2 mass flow", true);

      statistics.add_value("Bottom co2 mass flow",bottom_total_flux);
      statistics.set_precision ("Bottom co2 mass flow", 7);
      statistics.set_scientific ("Bottom co2 mass flow", true);

      statistics.add_value("Top co2 mass flow",top_total_flux);
      statistics.set_precision ("Top co2 mass flow", 7);
      statistics.set_scientific ("Top co2 mass flow", true);




      statistics.add_value("Left co2 solid flow",left_solid_flux);
      statistics.set_precision ("Left co2 solid flow", 7);
      statistics.set_scientific ("Left co2 solid flow", true);

      statistics.add_value("Right co2 solid flow",right_solid_flux);
      statistics.set_precision ("Right co2 solid flow", 7);
      statistics.set_scientific ("Right co2 solid flow", true);

      statistics.add_value("Bottom co2 solid flow",bottom_solid_flux);
      statistics.set_precision ("Bottom co2 solid flow", 7);
      statistics.set_scientific ("Bottom co2 solid flow", true);

      statistics.add_value("Top co2 solid flow",top_solid_flux);
      statistics.set_precision ("Top co2 solid flow", 7);
      statistics.set_scientific ("Top co2 solid flow", true);




      statistics.add_value("Left co2 fluid flow",left_fluid_flux);
      statistics.set_precision ("Left co2 fluid flow", 7);
      statistics.set_scientific ("Left co2 fluid flow", true);

      statistics.add_value("Right co2 fluid flow",right_fluid_flux);
      statistics.set_precision ("Right co2 fluid flow", 7);
      statistics.set_scientific ("Right co2 fluid flow", true);

      statistics.add_value("Bottom co2 fluid flow",bottom_fluid_flux);
      statistics.set_precision ("Bottom co2 fluid flow", 7);
      statistics.set_scientific ("Bottom co2 fluid flow", true);

      statistics.add_value("Top co2 fluid flow",top_fluid_flux);
      statistics.set_precision ("Top co2 fluid flow", 7);
      statistics.set_scientific ("Top co2 fluid flow", true);

      return std::pair<std::string, std::string> ("Mass fluxes through boundary parts:",
                                                  screen_text.str());
    }

    template <int dim>
    template <class Archive>
    void DetailedFluxStatistics<dim>::serialize (Archive &ar, const unsigned int)
    {
      ar &top_total_flux
      &bottom_total_flux
      &left_total_flux
      &right_total_flux
      &top_fluid_flux
      &bottom_fluid_flux
      &left_fluid_flux
      &right_fluid_flux
      &top_solid_flux
      &bottom_solid_flux
      &left_solid_flux
      &right_solid_flux;
    }


    template <int dim>
    void
    DetailedFluxStatistics<dim>::save (std::map<std::string, std::string> &status_strings) const
    {
      std::ostringstream os;

      // Serialize into a stringstream. Put the following into a code
      // block of its own to ensure the destruction of the 'oa'
      // archive triggers a flush() on the stringstream so we can
      // query the completed string below.
      {
        aspect::oarchive oa (os);
        oa << (*this);
      }

      status_strings["co2det"] = os.str();
    }


    template <int dim>
    void
    DetailedFluxStatistics<dim>::load (const std::map<std::string, std::string> &status_strings)
    {
      // see if something was saved
      if (status_strings.find("co2det") != status_strings.end())
        {
          std::istringstream is (status_strings.find("co2det")->second);
          aspect::iarchive ia (is);
          ia >> (*this);
        }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(DetailedFluxStatistics,
                                  "detailed flux statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the mass flux across boundaries. For each boundary "
                                  "indicator (see your geometry description for which boundary "
                                  "indicators are used), the mass flux is computed in outward "
                                  "direction, i.e., from the domain to the outside, using the "
                                  "formula $\\int_{\\Gamma_i} \\rho \\mathbf v \\cdot \\mathbf n$ "
                                  "where $\\Gamma_i$ is the part of the boundary with indicator $i$, "
                                  "$\\rho$ is the density as reported by the material model, "
                                  "$\\mathbf v$ is the velocity, and $\\mathbf n$ is the outward normal. "
                                  "\n\n"
                                  "As stated, this postprocessor computes the \\textit{outbound} mass "
                                  "flux. If you "
                                  "are interested in the opposite direction, for example from "
                                  "the core into the mantle when the domain describes the "
                                  "mantle, then you need to multiply the result by -1."
                                  "\n\n"
                                  ":::{note}\n"
                                  "In geodynamics, the term ``mass flux'' is often understood "
                                  "to be the quantity $\\rho \\mathbf v$, which is really a mass "
                                  "flux \\textit{density}, i.e., a vector-valued field. In contrast "
                                  "to this, the current postprocessor only computes the integrated "
                                  "flux over each part of the boundary. Consequently, the units of "
                                  "the quantity computed here are $\\frac{kg}{s}$.\n"
                                  ":::")
  }
}
