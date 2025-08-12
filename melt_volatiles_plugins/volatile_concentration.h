/*
  Copyright (C) 2013 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_volatiles_concentration_h
#define _aspect_postprocess_visualization_volatiles_concentration_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessor that takes an output vector
       * and computes a variable that represents the melt fraction at every
       * point.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class VolatilesConcentration
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          VolatilesConcentration ();

          std::vector<std::string>
          get_names () const override;

          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;


          UpdateFlags
          get_needed_update_flags () const override;

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

           /**
           * Read the parameters this class declares from the parameter file.
           */
          //double
          //calculate_liquid_solidus (double C, double Tm, double K, int solid);

        private:
          /**
           * Parameters for anhydrous melting of peridotite after Katz, 2003
           */
          double fluid_density_difference;       // Coefficients for T-dependence of distribution coefficients K^i
          double  rho_s;       // Coefficients for T-dependence of distribution coefficients K^i
          unsigned int porosity_idx;
          unsigned int ccl_idx;
          unsigned int ccs_idx;  
          unsigned int hcl_idx;
          unsigned int hcs_idx;  
          unsigned int mcl_idx;
          unsigned int mcs_idx;  
          unsigned int n_components;       // Coefficients for T-dependence of distribution coefficients K^i
          double cwt;
          double hwt;

      };
    }
  }
}

#endif
