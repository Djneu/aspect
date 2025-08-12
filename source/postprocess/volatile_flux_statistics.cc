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


#include <aspect/postprocess/volatile_flux_statistics.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/mesh_deformation/free_surface.h>
#include <aspect/melt.h>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    VolatileFluxStatistics<dim>::execute (TableHandler &statistics)
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

      std::map<types::boundary_id, double> local_co2_fluxes;
      std::map<types::boundary_id, double> local_h2o_fluxes;
      std::map<types::boundary_id, double> local_morb_fluxes;
      std::map<types::boundary_id, double> local_top_co2_fluxes;

      MaterialModel::MaterialModelInputs<dim> in(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(fe_face_values.n_quadrature_points, this->n_compositional_fields());
      
      MeltHandler<dim>::create_material_model_outputs(out);

      in.requested_properties = MaterialModel::MaterialProperties::density;
      MaterialModel::MeltOutputs<dim> *fluid_out = out.template get_additional_output<MaterialModel::MeltOutputs<dim>>();

      // Get the compositional values, and limit between 0 and 1.
      //double rho_s = 3200;
      //const double rho_l = 2700;

       // Get the boundary indicators of those boundaries with
      // a free surface.
      std::set<types::boundary_id> is_free_surface;
      if (this->get_parameters().mesh_deformation_enabled == true)
        is_free_surface = this->get_mesh_deformation_handler().get_free_surface_boundary_indicators();
      

      // for every surface face on which it makes sense to compute a
      // mass flux and that is owned by this processor,
      // integrate the normal mass flux given by the formula
      //   j =  \rho * v * n
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          for (const unsigned int f : cell->face_indices())
            if (cell->at_boundary(f))
              {

                const types::boundary_id boundary_ind
                  = cell->face(f)->boundary_id();
                  
                fe_face_values.reinit (cell, f);

                // Set use_strain_rates to false since we don't need viscosity
                in.reinit(fe_face_values, cell, this->introspection(), this->get_solution());

                this->get_material_model().evaluate(in, out);

                const FEValuesExtractors::Vector ex_u_f = this->introspection().variable("fluid velocity").extractor_vector();
                fe_face_values[ex_u_f].get_function_values (this->get_solution(), fluid_velocity_values);

                double local_co2_flux = 0;
                for (unsigned int q=0; q<fe_face_values.n_quadrature_points; ++q)
                  {
                  double cmorb_cl =  std::max(0.0, std::min(in.composition[q][ccl_idx],1.0));
                  double cmorb_cs =  std::max(0.0, std::min(in.composition[q][ccs_idx],1.0));
                  double Fvol =  std::max(0.0, std::min(in.composition[q][porosity_idx],1.0));

                  double rho_s = out.densities[q];
                  double rho_l = fluid_out->fluid_densities[q];

                  double avg_rho = Fvol*rho_l + (1 - Fvol)*rho_s;
                  double Fmass = Fvol*rho_l/avg_rho; 

                  local_co2_flux += (20./100 *  
                    (
                    (Fmass * cmorb_cl * rho_l * (fluid_velocity_values[q] * fe_face_values.normal_vector(q)))
                    + ((1 - Fmass) * cmorb_cs * rho_s * (in.velocity[q] * fe_face_values.normal_vector(q)))
                    )
                    * fe_face_values.JxW(q));     

              
                  // Old method using only solid velocity.
                  // (Fmass*cmorb_cl*rho_l + (1-Fmass)*cmorb_cs*rho_s) * 20/100
                  // * (in.velocity[q] * fe_face_values.normal_vector(q))
                  // * fe_face_values.JxW(q);

                  }                           

                const types::boundary_id boundary_indicator
                  = cell->face(f)->boundary_id();
                local_co2_fluxes[boundary_indicator] += local_co2_flux * in_years;                
              }

      

      // now communicate to get the global values
      std::map<types::boundary_id, double> global_co2_boundary_fluxes;
      {
        // first collect local values in the same order in which they are listed
        // in the set of boundary indicators
        const std::set<types::boundary_id>
        boundary_indicators
          = this->get_geometry_model().get_used_boundary_indicators ();

        std::vector<double> local_co2_values;
        std::vector<double> local_top_co2_values;
        local_co2_values.reserve(boundary_indicators.size());
        for (const auto p : boundary_indicators)
          local_co2_values.push_back (local_co2_fluxes[p]);

        // then collect contributions from all processors
        std::vector<double> global_co2_values (local_co2_values.size());
        Utilities::MPI::sum (local_co2_values, this->get_mpi_communicator(), global_co2_values);

        // and now take them apart into the global map again
        unsigned int index = 0;
        for (std::set<types::boundary_id>::const_iterator
             p = boundary_indicators.begin();
             p != boundary_indicators.end(); ++p, ++index)
              global_co2_boundary_fluxes[*p] = global_co2_values[index];

      }

      // now add all of the computed mass fluxes to the statistics object
      // and create a single string that can be output to the screen
      std::ostringstream screen_text;
      unsigned int index = 0;
      double top_co2_flux = 0;
      for (std::map<types::boundary_id, double>::const_iterator
           p = global_co2_boundary_fluxes.begin();
           p != global_co2_boundary_fluxes.end(); ++p, ++index)
           {

             // If it is a free surface we don't calculate outward flux.
            if (this->get_parameters().mesh_deformation_enabled == true)
              if(is_free_surface.find(p->first) != is_free_surface.end()) 
                 time_integrated_mass_flux += 0;
              else
                time_integrated_mass_flux += p->second * this->get_timestep() / year_in_seconds;
            else
              time_integrated_mass_flux += p->second * this->get_timestep() / year_in_seconds;

            // Find flux out top boundary, and convert back from ppm to kg/yr.
            // Do this regardless of whether it is a free surface.
            if(p->first == this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top"))
            {
              total_co2_degass += p->second * this->get_timestep() / year_in_seconds;
              top_co2_flux = p->second;
            }
           }


      // Now we find the global volatile mass

     // create a quadrature formula based on the compositional element alone.
      const Quadrature<dim> &quadrature_formula1 = this->introspection().quadratures.compositional_field_max;
      const unsigned int n_q_points = quadrature_formula1.size();

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula1,
                               update_values   |
                               update_quadrature_points |
                               update_gradients |
                               update_JxW_values);

      std::vector<double> Fvol_values(n_q_points);
      std::vector<double> ccl_values(n_q_points);
      std::vector<double> ccs_values(n_q_points);
      double local_co2_compositional_integrals = 0.0;

      MaterialModel::MaterialModelInputs<dim> in_fe(fe_values.n_quadrature_points, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out_fe(fe_values.n_quadrature_points, this->n_compositional_fields());
      
      MeltHandler<dim>::create_material_model_outputs(out_fe);

      in_fe.requested_properties = MaterialModel::MaterialProperties::density;
      MaterialModel::MeltOutputs<dim> *fluid_out_fe = out_fe.template get_additional_output<MaterialModel::MeltOutputs<dim>>();

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);

            // Set use_strain_rates to false since we don't need viscosity
            in_fe.reinit(fe_values, cell, this->introspection(), this->get_solution());

            this->get_material_model().evaluate(in_fe, out_fe);            

            fe_values[this->introspection().extractors.compositional_fields[porosity_idx]].get_function_values (this->get_solution(),
                Fvol_values);
            fe_values[this->introspection().extractors.compositional_fields[ccl_idx]].get_function_values (this->get_solution(),
                ccl_values);
            fe_values[this->introspection().extractors.compositional_fields[ccs_idx]].get_function_values (this->get_solution(),
                ccs_values);
                
            for (unsigned int q=0; q<n_q_points; ++q)
            {
              double rho_s = out_fe.densities[q];
              double rho_l = fluid_out_fe->fluid_densities[q];

              double avg_rho = Fvol_values[q]*rho_l + (1 - Fvol_values[q])*rho_s;
              double Fmass = Fvol_values[q]*rho_l/avg_rho;

              local_co2_compositional_integrals += (Fmass*ccl_values[q]*rho_l 
                                                  + (1-Fmass)*ccs_values[q]*rho_s ) * 20/100
                                                  *fe_values.JxW(q);
            }


          }

      // compute the sum over all processors
      const double global_co2_compositional_integrals =
      Utilities::MPI::sum (local_co2_compositional_integrals,
                           this->get_mpi_communicator());
      
      // Positive indicates outward flow, so to see the total amount of mass we have had,
      // we add that to the value.
      double conserve_co2 = global_co2_compositional_integrals + time_integrated_mass_flux;

      // Co2
      statistics.add_value("Total co2 mass flow",-time_integrated_mass_flux);
      statistics.set_precision ("Total co2 mass flow", 7);
      statistics.set_scientific ("Total co2 mass flow", true);

      statistics.add_value("Global co2 mass",global_co2_compositional_integrals);
      statistics.set_precision ("Global co2 mass", 7);
      statistics.set_scientific ("Global co2 mass", true);

      statistics.add_value("Total co2",conserve_co2);
      statistics.set_precision ("Total co2", 7);
      statistics.set_scientific ("Total co2", true);

      statistics.add_value("Co2 degass rate",top_co2_flux);
      statistics.set_precision ("Co2 degass rate", 7);
      statistics.set_scientific ("Co2 degass rate", true);

      statistics.add_value("Total Co2 degass",total_co2_degass);
      statistics.set_precision ("Total Co2 degass", 7);
      statistics.set_scientific ("Total Co2 degass", true);
      

      return std::pair<std::string, std::string> ("Writing volatile statistics",
                                                screen_text.str());

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(VolatileFluxStatistics,
                                  "volatile statistics",
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
