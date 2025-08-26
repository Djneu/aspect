
/*
  Copyright (C) 2019 - 2021 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_pore_fluid_pressure_h
#define _aspect_material_model_pore_fluid_pressure_h

#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>

#include<deal.II/fe/component_mask.h>
#include <deal.II/matrix_free/fe_point_evaluation.h>

#include <aspect/material_model/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace Rheology
    {
      template <int dim>
      class PoreFluidPressure : public ::aspect::SimulatorAccess<dim>
      {
        public:
          /**
           * Declare the parameters this function takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm,
                            const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition = nullptr);

          /**
           * A function that computes by how much the rheologic parameters change
           * if strain weakening is applied. Given a compositional field with
           * the index j and a vector of all compositional fields, it returns
           * reduction factors for the cohesion, friction angle and the prefactor
           * of the viscous flow law(s) used in the computation for that composition.
           */
          double
          compute_fluid_ratio(const unsigned int composition,
                              const double depth) const;

        private:

          std::vector<double> fluid_ratio_top;
          std::vector<double> fluid_ratio_base;
          std::vector<double> fluid_ratio_length;
          double fluid_ratio_cutoff;
          double fluid_cutoff_taper;
          bool use_pore_fluid_pressure;
      };
    }
  }
}
#endif
