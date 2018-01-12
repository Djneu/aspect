/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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

#include <aspect/postprocess/visualization/phase_tracker.h>
#include <aspect/material_model/phase_transitions.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      template <int dim>
      PhaseTracker<dim>::
      PhaseTracker ()
        :
        DataPostprocessorScalar<dim> ("phase_tracker",
                                      update_values | update_q_points )
      {}



      template <int dim>
      void
      PhaseTracker<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {

        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        // in case the material model computes the melt fraction iself, we use that output
        if (const MaterialModel::PhaseTransitions<dim> *
            melt_material_model = dynamic_cast <const MaterialModel::PhaseTransitions<dim>*> (&this->get_material_model()))
          {
            MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                       this->introspection());
            MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                         this->n_compositional_fields());

            // Compute the melt fraction...
            this->get_material_model().evaluate(in, out);

            std::vector<double> phase_tracker(n_quadrature_points);
            melt_material_model->phase_tracker(in, phase_tracker);

            for (unsigned int q=0; q<n_quadrature_points; ++q)
              computed_quantities[q](0) = phase_tracker[q];
          }

       /*onst unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        Assert (computed_quantities[0].size() == 1,                   ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,           ExcInternalError());

        // Set use_strain_rates to false since we have no need for viscosity.
        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection(), false);
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        std::vector<double> phase_tracker(n_quadrature_points);
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          computed_quantities[q](0) = phase_tracker[q];*/
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
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(PhaseTracker,
                                                  "phase tracker",
                                                  "A visualization output object that generates output "
                                                  "for the thermal conductivity $k$.")
    }
  }
}
