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

#ifndef _aspect_material_reaction_co2_melt_h
#define _aspect_material_reaction_co2_melt_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    namespace ReactionModel
    {

      /**
      * A melt model that calculates melt fraction and entropy change
      * according to the melting model for dry peridotite of Katz, 2003.
      * This also includes a computation of the latent heat of melting (if the latent heat
      * heating model is active).
      *
      * These functions can be used in the calculation of melting and melt transport
      * in the melt_simple material model and can be extended to other material models
      *
      * @ingroup ReactionModel
      */
      template <int dim>
      class Co2Melt : public ::aspect::SimulatorAccess<dim>
      {
        public:
          // constructor
          Co2Melt();

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
          parse_parameters (ParameterHandler &prm);

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
                              bool compute_solidus,
                              const double temp) const;


          /**
           * Compute all the reaction rate variables needed for a reactive transport model based on the
           * Katz 2003 formulation. Takes the material model inputs @p in to compute the material model outputs @p out.
           * This function mainly fills the reaction_rate_out object but populates out.reaction_terms,
           * out.entropy_derivative_pressure and entropy_derivative_temperature
          */
          void calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                               typename Interface<dim>::MaterialModelOutputs &out) const;

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
          double rho_l;       // Coefficients for T-dependence of distribution coefficients K^i
          double  rho_s;       // Coefficients for T-dependence of distribution coefficients K^i
          double melting_time_scale;
      };
    }

  }
}

#endif
