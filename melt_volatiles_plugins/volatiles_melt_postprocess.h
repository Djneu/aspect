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


#ifndef _aspect_postprocess_visualization_volatiles_melt_post_h
#define _aspect_postprocess_visualization_volatiles_melt_post_h

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
      class VolatilesMeltPost
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          VolatilesMeltPost ();

          std::vector<std::string>
          get_names () const override;

          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const override;


          UpdateFlags
          get_needed_update_flags () const override;

          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;

          // Calculate melting temperature for each composition.
          std::vector<double>
          melting_temperatures(const double pressure) const;

          // Calculate melting temperature for each composition.
          double
          compute_residual(std::vector<double> composition,
                      std::vector<double> K,
                      bool compute_solidus) const;

          // Calculate melting temperature for each composition.
          std::vector<double>
          partition_coefficients(const double pressure,
                               const double temperature) const;

          double
          T_solidus_liquidus (const double pressure, 
                              std::vector<double> composition, 
                              bool compute_solidus) const;

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

          std::vector<double> T0;      // Pure component melting points at P=0
          std::vector<double> A;       // Coefficients for linear P-dependence of T_m^i
          std::vector<double> B;       // Coefficients for quadratic P-dependence of T_m^i
          std::vector<double> L;       // Latent heat of pure components
          std::string K_T_mode;        // Type of parameterization for K^i(T)
          std::vector<double> R;       // Coefficients for T-dependence of distribution coefficients K^i
          unsigned int n_components;       // Coefficients for T-dependence of distribution coefficients K^i
          double porosity;       // Coefficients for T-dependence of distribution coefficients K^i
          double reference_darcy_coefficient () const;
          double rho_l;       // Coefficients for T-dependence of distribution coefficients K^i
          double  rho_s;       // Coefficients for T-dependence of distribution coefficients K^i

          unsigned int melt_idx;
          unsigned int mcl_idx;
          unsigned int mcs_idx;
          unsigned int ccl_idx;
          unsigned int ccs_idx;  
          unsigned int hcl_idx;
          unsigned int hcs_idx;  
          double pressure_max;
      };
    }
  }
}

#endif
