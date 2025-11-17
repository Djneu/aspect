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

#ifndef _aspect_material_reaction_volatiles_melt_h
#define _aspect_material_reaction_volatiles_melt_h

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
      class VolatilesMelt : public ::aspect::SimulatorAccess<dim>
      {
        public:
          // constructor
          VolatilesMelt();

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
          small_vector<double>
          melting_temperatures(const double pressure,
                              const std::vector<double> &A,
                              const std::vector<double> &B,
                              const std::vector<double> &L,
                              const std::vector<double> &T0) const;

          // Calculate melting temperature for each composition.
          double
          compute_residual(const std::vector<double> &composition,
                      const small_vector<double> &K,
                      const bool compute_solidus) const;

          // Calculate melting temperature for each composition.
          small_vector<double>
          partition_coefficients(const double pressure,
                               const double temperature,
                               const std::vector<double> &A,
                               const std::vector<double> &B,
                               const std::vector<double> &L,
                               const std::vector<double> &T0,
                               const std::vector<double> &R) const;

          double
          T_solidus_liquidus (const double pressure, 
                              std::vector<double> composition, 
                              bool compute_solidus,                              
                              std::vector<double> A,
                              std::vector<double> B,
                              std::vector<double> L,
                              std::vector<double> T0,
                              std::vector<double> R) const;


          /**
           * Compute all the reaction rate variables needed for a reactive transport model based on the
           * Katz 2003 formulation. Takes the material model inputs @p in to compute the material model outputs @p out.
           * This function mainly fills the reaction_rate_out object but populates out.reaction_terms,
           * out.entropy_derivative_pressure and entropy_derivative_temperature
          */
          void calculate_reaction_rate_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                               typename Interface<dim>::MaterialModelOutputs &out) const;

                              

          /**
           * Compute all the fluid variables needed for a reactive transport model based on the
           * Katz 2003 formulation. This function fills melt outputs, the out object should already contain
           * outputs for the solid and this function uses the inputs @p in and the solid outputs @p out
           * to fill MeltOutputs. Solid outputs such as out.Thermal_expansion_coefficients are expected
           * to have already been computed when this function is called. Solid viscosities are also modified
           * in the out object here because the presence of melt weakens the material.
           */
          void calculate_fluid_outputs(const typename Interface<dim>::MaterialModelInputs &in,
                                       typename Interface<dim>::MaterialModelOutputs &out,
                                       const double reference_T) const;


        /**
         * Function that follows Keller and Katz, 2016, to calculate equilibrium
         * solid and liquid components, and melt fraction for a multi-component system.
         * This function returns the volume melt fraction, melt reaction rate,
         * and solid and liquid reaction rates for each component.
         */
        std::tuple<double, double, std::vector<double>, std::vector<double>, double>
        equilibrium (std::vector<double> composition, 
                     const double temperature, 
                     const double pressure,
                     const double depth,
                     const double rho_ss,
                     const int ep,
                     const double x) const;

          double reference_darcy_coefficient () const;

        /**
         * Make sure the reaction rate is corrected that any value already outside of 
         * 0 and 1, or outside of 0 and 1 with the reaction change, does not leave bounds.
         */
        virtual
        double
        get_reaction_rate (const double old_value,
                           const double reaction_rate,
                           const double time_step,
                           const double depth,
                           const double x) const;

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          initialize_component_values (std::vector<double> &Ac,
                                                          std::vector<double> &Bc,
                                                          std::vector<double> &Lc,
                                                          std::vector<double> &Tc,
                                                          std::vector<double> &Rc,
                                                          const double pressure) const;

          double 
          linear_interpolation(double P, double P1, double P2, double v1, double v2) const;

          

        private:
          /**
          * Parameters for reaction model
          */

          std::vector<double> T0i;      // Pure component melting points at P=0
          std::vector<double> Ai;       // Coefficients for linear P-dependence of T_m^i
          std::vector<double> Bi;       // Coefficients for quadratic P-dependence of T_m^i
          std::vector<double> Li;       // Latent heat of pure components
          std::string K_T_mode;        // Type of parameterization for K^i(T)
          std::vector<double> Ri;       // Coefficients for T-dependence of distribution coefficients K^i
          unsigned int n_components;       // Coefficients for T-dependence of distribution coefficients K^i
          double pressure_max;
          double fluid_density_difference;       // Coefficients for T-dependence of distribution coefficients K^i

          unsigned int porosity_idx;  // porosity index
          unsigned int mcl_idx;       // liquid morb index
          unsigned int mcs_idx;       // solid morb index
          unsigned int ccl_idx;       // liquid carbonated morb
          unsigned int ccs_idx;       // solid carbonated morb
          unsigned int hcl_idx;       // liquid hydrated morb
          unsigned int hcs_idx;       // solid hydrated morb          

          double xi_0;                               // rock viscosity constant
          double viscosity_fluid;                    // melt viscosity constant
          double thermal_bulk_viscosity_exponent;    
          double alpha_phi;                          // melt weakening factor
          double melt_compressibility;               
          double melting_time_scale;                 // reaction time
          double melt_bulk_modulus_derivative;       //
          double reference_permeability;             // permeability constant 
          double maximum_permeability;             // permeability constant   
          double compaction_viscosity_ratio;
          bool use_fractional_melting;
          double extraction_depth;
          double extraction_width;
          double use_extraction_patch;
          bool use_simons_law;

          mutable int totaln = 0;
          mutable double totaltime = 0;
      };
    }

  }
}

#endif
