/*
Copyright (C) 2011 - 2017 by the authors of the ASPECT code.

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


#include <aspect/material_model/phase_transitions.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>

namespace aspect
{
  namespace MaterialModel
  {
    template <int dim>
    void
    PhaseTransitions<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.position.size(); ++i)
      {
        const double temperature = in.temperature[i];
        const double pressure = in.pressure[i];
        const Point<dim> position = in.position[i];
        double depth = this->get_geometry_model().depth(position);

        //constant properties
        out.compressibilities[i] = reference_compressibility;
        out.specific_heat[i] = reference_specific_heat;

        /*thermal conductivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
         For now only uses perovskite periclase phase to calculate.*/
        out.thermal_conductivities[i] = (d0+(d1*depth))*pow((300/temperature),d2);

        /*thermal expansivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
         For now only uses fperovskite periclase phase to calculate.*/
        out.thermal_expansion_coefficients[i] = (b0+(b1*temperature)+b2/(temperature*temperature))*exp(-b3*depth);

        //placeholder viscosity equation
        if (temperature !=0.0)
          {
            out.viscosities[i] =B/2*exp(activation_energy/(8.314*temperature));
          }



        //density equation with pressure and temperature dependence, likely will change when adiabatic conditions are introduced.
        double density_temperature_dependence = 1.0;
        if (this->include_adiabatic_heating ())
        {
          density_temperature_dependence -= thermal_alpha * (temperature - this->get_adiabatic_conditions().temperature(position));
        }
        else
          density_temperature_dependence -= thermal_alpha * (temperature-reference_T);

        const double density_pressure_dependence = reference_rho * out.compressibilities[i] * (pressure - this->get_surface_pressure());
        out.densities[i] = (reference_rho + density_pressure_dependence)*density_temperature_dependence;
      }
    }


    template <int dim>
    double
    PhaseTransitions<dim>::
    reference_viscosity () const
    {
      return eta;
    }

    template <int dim>
    bool
    PhaseTransitions<dim>::
    is_compressible () const
    {
      if (reference_compressibility > 0)
        return true;
      else
        return false;
    }



    template <int dim>
    void
    PhaseTransitions<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Phase transitions");
          {
            prm.declare_entry ("Reference density", "3300",
                               Patterns::Double (0),
                               "Reference density $\\rho_0$. Units: $kg/m^3$.");
            prm.declare_entry ("Reference temperature", "293",
                               Patterns::Double (0),
                               "The reference temperature $T_0$. The reference temperature is used "
                               "in the density formula. Units: $K$.");
            prm.declare_entry ("Viscosity", "5e24",
                               Patterns::Double (0),
                               "The value of the viscosity $\\eta$. Units: $kg/m/s$.");
            prm.declare_entry ("Thermal conductivity", "4.7",
                               Patterns::Double (0),
                               "The value of the thermal conductivity $k$. "
                               "Units: $W/m/K$.");
            prm.declare_entry ("Reference specific heat", "1250",
                               Patterns::Double (0),
                               "The value of the specific heat $C_p$. "
                               "Units: $J/kg/K$.");
            prm.declare_entry ("Thermal expansion coefficient", "4e-5",
                               Patterns::Double (0),
                               "The value of the thermal expansion coefficient $\\beta$. "
                               "Units: $1/K$.");
            prm.declare_entry ("Compressibility", "5.124e-12",
                               Patterns::Double (0),
                               "The value of the compressibility $\\kappa$. "
                               "Units: $1/Pa$.");
            prm.declare_entry ("Viscosity constant", "9.1e11",
                               Patterns::Double(0),
                               "Preexponential constant for viscosity." "Units: None");
            prm.declare_entry ("Activation energy", "300000",
                               Patterns::Double(0),
                               "Activation energy for viscosity equation." "Units: J/mol");
            prm.declare_entry ("b0", "2.68e-5",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/K");
            prm.declare_entry ("b1", "2.77e-9",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/K^2");
            prm.declare_entry ("b2", "-1.21",
                               Patterns::Double(),
                               "coefficient for depth dependent thermal expansivity" "Units: K"); 
            prm.declare_entry ("b3", "3.76e-7",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/m");
            prm.declare_entry ("d0", "3.48",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal conductivity" "Units: Wm^-1K^-1");
            prm.declare_entry ("d1", "5.17e-6",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal conductivity" "Units: Wm^-2K^-1");
            prm.declare_entry ("d2", "0.31",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal conductivity" "Units: None");

          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }

    template <int dim>
    void
    PhaseTransitions<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Phase transitions");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          reference_compressibility  = prm.get_double ("Compressibility");
          eta                        = prm.get_double ("Viscosity");
          k_value                    = prm.get_double ("Thermal conductivity");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          B                          = prm.get_double ("Viscosity constant");
          activation_energy          = prm.get_double ("Activation energy");
          b0                         = prm.get_double ("b0");
          b1                         = prm.get_double ("b1");
          b2                         = prm.get_double ("b2");
          b3                         = prm.get_double ("b3");
          d0                         = prm.get_double ("d0");
          d1                         = prm.get_double ("d1");
          d2                         = prm.get_double ("d2");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::none;
      this->model_dependence.density = NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::none;
      this->model_dependence.thermal_conductivity = NonlinearDependence::none;
    }
  }
}

namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(PhaseTransitions,
                                   "Phase transitions",
                                   "Material model to eventually have phase transitions and melt.")
  }
}

