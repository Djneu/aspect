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
    //temperature and pressure dependent viscosity, set up similar to damage_rheology material model
    template <int dim>
    double
    PhaseTransitions<dim>::
    calculate_viscosity ( const double &pressure,
                          const double &temperature,
                          const Point<dim> &position,
                          const SymmetricTensor<2,dim> &strain_rate) const
    {
      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      double current_viscosity = eta;
      double vdis              = eta;

      const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                        ?
                                        this->get_adiabatic_conditions().pressure(position)
                                        :
                                        pressure;

      //energy term for viscosity calculated with diffusion
      double energy_term = exp((activation_energy[0] + activation_volume[0] * adiabatic_pressure)
                         / (1 * constants::gas_constant * temperature));

      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = exp((activation_energy[0] + activation_volume[0] * adiabatic_pressure)
              / (1 * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > 1e6)
            energy_term = adiabatic_energy_term * 1e6;
          if (temperature_dependence < 1.0 / 1e6)
            energy_term = adiabatic_energy_term / 1e6;
        }

      const double strain_rate_dependence = (1.0 - 1) / 1;

      double vdiff = pow(A[0],-1.0/1)
             * std::pow(second_strain_rate_invariant,strain_rate_dependence)
             * pow(grain_size, 3/1)
             * energy_term;


      //energy term for calculation of viscosity from dislocation creep
      double energy_term_dis = exp((activation_energy[1] + activation_volume[1] * adiabatic_pressure)
                         / (3.5 * constants::gas_constant * temperature));
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term_dis
            = exp((activation_energy[1] + activation_volume[1] * adiabatic_pressure)
              / (3.5 * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double temperature_dependence_dis = energy_term_dis / adiabatic_energy_term_dis;
          if (temperature_dependence_dis > 1e6)
            energy_term_dis = adiabatic_energy_term_dis * 1e6;
          if (temperature_dependence_dis < 1.0 / 1e6)
            energy_term_dis = adiabatic_energy_term_dis / 1e6;
        }

      const double strain_rate_dependence_dis = (1.0 - 3.5) / 3.5;

      if(second_strain_rate_invariant !=0.0)
      {
             vdis = pow(A[1],-1.0/3.5)
             * std::pow(second_strain_rate_invariant,strain_rate_dependence_dis)
             * energy_term_dis;
      }
      else
      {
             vdis = pow(A[1],-1.0/3.5)
             * std::pow(1e-15,strain_rate_dependence_dis)
             * energy_term_dis;
      }

      if(std::abs(second_strain_rate_invariant) > 1e-30)
        current_viscosity = vdis * vdiff / (vdis + vdiff);
      else
        current_viscosity = vdiff;

     return current_viscosity;
    }


    template <int dim>
    void
    PhaseTransitions<dim>::
    evaluate(const MaterialModel::MaterialModelInputs<dim> &in,
             MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      for (unsigned int i=0; i < in.temperature.size(); ++i)
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


        if(in.strain_rate.size())
        {
                      out.viscosities[i] = std::min(std::max(min_eta,calculate_viscosity(pressure,
                                                                       temperature,
                                                                       position,
                                                                       in.strain_rate[i])),max_eta);
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

            //viscosity parameters
            prm.declare_entry ("Viscosity prefactor", "9.1e11",
                               Patterns::List (Patterns::Double(0)),
                               "Viscosity prefactor with first entry being diffusion, second dislocation" "Units: None");
            prm.declare_entry ("Activation energy", "300000",
                               Patterns::List (Patterns::Double(0)),
                               "Activation energy for viscosity equation." "Units: J/mol");
            prm.declare_entry ("Activation volume", "6e-6",
                               Patterns::List (Patterns::Double(0)),
                               "Activation volume for viscosity equation." "Units: m^3/mol");
            prm.declare_entry ("Maximum viscosity","1.0e24",Patterns::Double(0),
                               "Reference strain rate for first time step. Units: Kg/m/s");
            prm.declare_entry ("Minimum viscosity", "1.0e18", Patterns::Double(0),
                               "Stabilizes strain dependent viscosity. Units: Kg/m/s");
            prm.declare_entry ("Grain size", "0.003", Patterns::Double(0),
                               "Stabilizes strain dependent viscosity. Units: m");

            //thermal expansivity and conductivity parameters
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

          //viscosity parameters
          max_eta                    = prm.get_double ("Maximum viscosity");
          min_eta                    = prm.get_double("Minimum viscosity");
          grain_size                 = prm.get_double("Grain size");
          A                          = Utilities::string_to_double
                                       (Utilities::split_string_list(prm.get ("Viscosity prefactor")));
          activation_energy          = Utilities::string_to_double
                                       (Utilities::split_string_list(prm.get ("Activation energy")));
          activation_volume          = Utilities::string_to_double
                                       (Utilities::split_string_list(prm.get ("Activation volume")));


          if(A.size() !=2 || activation_volume.size() !=2 || activation_energy.size() !=2)
               AssertThrow(false, ExcMessage("Error: Too many inputs for activation volume, activation energy, or viscosity prefactor."));

          //thermal exansivity and conductivity parameters
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
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure;
      this->model_dependence.density = NonlinearDependence::temperature;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::temperature;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature;
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

