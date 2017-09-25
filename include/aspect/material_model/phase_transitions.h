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

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /*Beginning of material model to test viscosity and density changes.
    */
    template <int dim>
    class PhaseTransitions : public MaterialModel::Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:

        virtual void evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
                              MaterialModel::MaterialModelOutputs<dim> &out) const;

        virtual bool is_compressible () const;
        virtual double reference_viscosity () const;

        static
        void
        declare_parameters (ParameterHandler &prm);

        virtual
        void
        parse_parameters (ParameterHandler &prm);

      private:
        double reference_rho;
        double reference_T;
        double reference_compressibility;
        double reference_specific_heat;

        //viscosity variables
        double eta;
        double max_eta;
        double min_eta;
        std::vector<double> grain_size;
        std::vector<double> diffusion_activation_energy;
        std::vector<double> diffusion_activation_volume;
        std::vector<double> diffusion_prefactor;
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
        std::vector<double> phase_prefactors;

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
        phase_function_derivative (const Point<dim> &position,
                                   const double temperature,
                                   const double pressure,
                                   const int phase) const;
        virtual
        unsigned int
        get_phase_index (const Point<dim> &position,
                         const double temperature,
                         const double pressure) const;
        
    };
  }
}

#endif
