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
        double thermal_alpha;

        //viscosity variables
        double eta;
        double max_eta;
        double min_eta;
        double grain_size;
        std::vector<double> activation_energy;
        std::vector<double> activation_volume;
        std::vector<double> A;

        //thermal expansivity and conductivity variables
        double k_value;
        double b0;
        double b1;
        double b2;
        double b3;
        double d0; 
        double d1; 
        double d2;

        virtual
        double
        calculate_viscosity ( const double &pressure,
                              const double &temperature,
                              const Point<dim> &position,
                              const SymmetricTensor<2,dim> &strain_rate) const;
        
    };
  }
}

#endif
