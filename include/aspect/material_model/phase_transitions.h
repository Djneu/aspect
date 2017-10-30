/*
  Copyright (C) 2013 - 2017 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_phase_transitions_h
#define _aspect_material_model_phase_transitions_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/melt.h>

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /*Beginning of material model to test viscosity and density changes.
    */
    template <int dim>
    class PhaseTransitions : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>, public MaterialModel::MeltFractionModel<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        virtual bool is_compressible () const;
        virtual double reference_viscosity () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        double test = 0.;

        virtual
        void
        parse_parameters (ParameterHandler &prm);

        /**
         * @}
         */

        virtual void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                                     std::vector<double> &melt_fractions) const;

      private:
        double reference_rho;
        double reference_T;
        double reference_compressibility;
        double reference_specific_heat;

 
        // Start of melting parameters

        // Melt thermal expansivity
        double melt_thermal_alpha;

        /**
         * Parameters for anhydrous melting of peridotite after Katz, 2003
         */

        // for the solidus temperature
        double A1;   // °C
        double A2; // °C/Pa
        double A3; // °C/(Pa^2)

        // for the lherzolite liquidus temperature
        double B1;   // °C
        double B2;   // °C/Pa
        double B3; // °C/(Pa^2)

        // for the liquidus temperature
        double C1;   // °C
        double C2;  // °C/Pa
        double C3; // °C/(Pa^2)

        // for the reaction coefficient of pyroxene
        double r1;     // cpx/melt
        double r2;     // cpx/melt/GPa
        double M_cpx;  // mass fraction of pyroxene

        // melt fraction exponent
        double beta;

        // entropy change upon melting
        double peridotite_melting_entropy_change;

        // the relative density of molten material (compared to solid)
        double relative_melt_density;

        /**
         * Percentage of material that is molten. Melting model after Katz,
         * 2003 (for peridotite) and Sobolev et al., 2011 (for pyroxenite)
         */
        virtual
        double
        melt_fraction (const double temperature,
                       const double pressure,
                       const std::vector<double> &compositional_fields,
                       const Point<dim> &position) const;

        virtual
        double
        peridotite_melt_fraction (const double temperature,
                                  const double pressure,
                                  const std::vector<double> &compositional_fields,
                                  const Point<dim> &position) const;
        double
        entropy_derivative_melt ( const double temperature,
                                  const double pressure,
                                  const std::vector<double> &compositional_fields,
                                  const Point<dim> &position,
                                  const NonlinearDependence::Dependence dependence) const;

        // grain evolution parameters
        double gas_constant; // J/(K*mol)
        std::vector<double> grain_growth_activation_energy;
        std::vector<double> grain_growth_activation_volume;
        std::vector<double> grain_growth_rate_constant;
        std::vector<double> grain_growth_exponent;
        std::vector<double> reciprocal_required_strain;
        std::vector<double> recrystallized_grain_size;
        double min_grain_size;
        double pv_grain_size_scaling;

        bool advect_log_gransize;

        // for paleowattmeter
        bool use_paleowattmeter;
        std::vector<double> grain_boundary_energy;
        std::vector<double> boundary_area_change_work_fraction;
        std::vector<double> geometric_constant;

        //viscosity variables
        double dislocation_viscosity_iteration_threshold;
        unsigned int dislocation_viscosity_iteration_number;
        std::vector<double> dislocation_creep_exponent;
        double max_temperature_dependence_of_eta;
        double eta;
        double max_eta;
        double min_eta;
        std::vector<double> constant_grain_size;
        std::vector<double> diffusion_activation_energy;
        std::vector<double> diffusion_activation_volume;
        std::vector<double> diffusion_prefactor;
        std::vector<double> diffusion_creep_grain_size_exponent;
        std::vector<double> dislocation_activation_energy;
        std::vector<double> dislocation_activation_volume;
        std::vector<double> dislocation_prefactor;
 
       //thermal expansivity and conductivity variables
        double k_value;
        double thermal_alpha;
        std::vector<double> a0;
        std::vector<double> a1;
        std::vector<double> a2;
        std::vector<double> a3;
        std::vector<double> c0; 
        std::vector<double> c1; 
        std::vector<double> c2;

        //phase boundary variables
        std::vector<double> transition_depths;
        std::vector<double> transition_temperatures;
        std::vector<double> transition_widths;
        std::vector<double> transition_slopes;
        std::vector<double> density_jumps;
        std::vector<double> temperature_jumps;  
        std::vector<double> phase_prefactors;

        /**
         * Function that takes an object in the same format
         * as in.composition as argument and converts the
         * vector that corresponds to the grain size to its
         * logarithms and back and limits the grain size to
         * a global minimum.
         * @in normal_to_log: if true, convert from the grain
         * size to its logarithm, otherwise from log to grain
         * size
         */
        virtual
        void
        convert_log_grain_size (const bool normal_to_log,
                                std::vector<double> &compositional_fields) const;

        /**
         * Rate of grain size growth (Ostwald ripening) or reduction
         * (due to phase transformations) in dependence on temperature
         * pressure, strain rate, mineral phase and creep regime.
         * We use the grain size evolution laws described in Solomatov
         * and Reese, 2008. Grain size variations in the Earth’s mantle
         * and the evolution of primordial chemical heterogeneities,
         * J. Geophys. Res., 113, B07408.
         */
        virtual
        double
        grain_size_growth_rate (const double                  temperature,
                                const double                  pressure,
                                const std::vector<double>    &compositional_fields,
                                const SymmetricTensor<2,dim> &strain_rate,
                                const Tensor<1,dim>          &velocity,
                                const Point<dim>             &position,
                                const unsigned int            phase_index,
                                const int                     crossed_transition) const;

        virtual double viscosity (const double                  temperature,
                                  const double                  pressure,
                                  const std::vector<double>    &compositional_fields,
                                  const SymmetricTensor<2,dim> &strain_rate,
                                  const Point<dim>             &position) const;


        virtual double diffusion_viscosity (const double      temperature,
                                            const double      pressure,
                                            const std::vector<double>    &compositional_fields,
                                            const SymmetricTensor<2,dim> &,
                                            const Point<dim> &position) const;

        /**
         * This function calculates the dislocation viscosity. For this purpose
         * we need the dislocation component of the strain rate, which we can
         * only compute by knowing the dislocation viscosity. Therefore, we
         * iteratively solve for the dislocation viscosity and update the
         * dislocation strain rate in each iteration using the new value
         * obtained for the dislocation viscosity. The iteration is started
         * with a dislocation viscosity calculated for the whole strain rate
         * unless a guess for the viscosity is provided, which can reduce the
         * number of iterations significantly.
         */
        virtual double dislocation_viscosity (const double      temperature,
                                              const double      pressure,
                                              const std::vector<double>    &compositional_fields,
                                              const SymmetricTensor<2,dim> &strain_rate,
                                              const Point<dim> &position,
                                              const double viscosity_guess = 0) const;

        /**
         * This function calculates the dislocation viscosity for a given
         * dislocation strain rate.
         */
        double dislocation_viscosity_fixed_strain_rate (const double      temperature,
                                                        const double      pressure,
                                                        const std::vector<double> &,
                                                        const SymmetricTensor<2,dim> &dislocation_strain_rate,
                                                        const Point<dim> &position) const;

        virtual
        double
        calculate_viscosity ( const double &pressure,
                              const double &temperature,
                              const Point<dim> &position,
                              const SymmetricTensor<2,dim> &strain_rate,
                              const int phase) const;

        virtual
        double
        phase_function (const Point<dim> &position,
                        const double temperature,
                        const double pressure,
                        const int phase) const;

        virtual
        double
        Pphase_function_derivative (const Point<dim> &position,
                                   const double temperature,
                                   const double pressure,
                                   unsigned int phase) const;

        virtual
        double
        Tphase_function_derivative (const Point<dim> &position,
                                   const double temperature,
                                   const double pressure,
                                   unsigned int phase) const;

        virtual
        unsigned int
        get_phase_index (const Point<dim> &position,
                         const double temperature,
                         const double pressure) const;
        
    };
  }
}

#endif
