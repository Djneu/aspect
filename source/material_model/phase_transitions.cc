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
      double vdiff=eta;

      return vdiff;
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    phase_function (const Point<dim> &,
                    const double temperature,
                    const double pressure,
                    const int phase) const
    {

          double transition_pressure;
          double pressure_width;
          double width_temp;
          double transition_temp;
          double phase_func;

          //converting depth and width given to pressures
          const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
          const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + 0.5 * transition_widths[phase]);
          const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - 0.5 * transition_widths[phase]);
          transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
          pressure_width = (this->get_adiabatic_conditions().pressure(transition_plus_width)
                             - this->get_adiabatic_conditions().pressure(transition_minus_width));
          width_temp = transition_widths[phase];



          //checking to use given transition temperature or find it from adiabatic temperature
          if(auto_temp)
          {
            transition_temp = this->get_adiabatic_conditions().temperature(transition_point);
          }
          else
            transition_temp = transition_temperatures[phase];

          //std::ofstream mfile;
	  //mfile.open("intpressure.txt", std::ios::app);
	 // mfile<<transition_depths[phase]<<"  "<<transition_pressure<<"  "<<transition_temp<<std::endl;
	 // mfile.close();


          // then calculate the deviation from the transition point (both in temperature
          // and in pressure)
          double pressure_deviation = pressure - transition_pressure
                                      - transition_slopes[phase] * (temperature - transition_temp);

          // last, calculate the percentage of material that has undergone the transition
          // (also in dependence of the phase transition width - this is an input parameter)
          // use delta function for width = 0
          if (width_temp==0)
            (pressure_deviation > 0) ? phase_func = 1 : phase_func = 0;
          else
            phase_func = 0.5*(1.0 + std::tanh(pressure_deviation / pressure_width));
          return phase_func;  
    }

    template <int dim>
    unsigned int
    PhaseTransitions<dim>::
    get_phase_index (const Point<dim> &position,
                     const double temperature,
                     const double pressure) const
    {
      unsigned int phase_index = 0;
      if(transition_depths.size()>0)
        if(phase_function(position, temperature, pressure, transition_depths.size()-1) == 1)
          phase_index = transition_depths.size();

      for(unsigned int j=1;j<transition_depths.size();++j)
        if(phase_function(position, temperature, pressure, j) != phase_function(position, temperature, pressure, j-1))
          phase_index = j;


          //std::ofstream mfile;
	  //mfile.open("stuff.txt", std::ios::app);
	// mfile<<depth<<"  "<<phase_index<<std::endl;
	  //mfile.close();


       return phase_index;
    }
         


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
        unsigned int number_of_phase_transitions = transition_depths.size();


        //constant properties
        out.compressibilities[i] = reference_compressibility;
        out.specific_heat[i] = reference_specific_heat;

        /*thermal conductivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
         For now only uses perovskite periclase phase to calculate.*/
        if(k_value == 0.0 && this->get_adiabatic_conditions().is_initialized())
        {
            const unsigned int ol_index = get_phase_index(position, temperature, pressure);
           // out.thermal_conductivities[i] = (c0[1]+(c1[1]*pressure*1e-9))*pow((300/temperature),c2[1]);
            out.thermal_conductivities[i] = c0[ol_index];
         //std::ofstream mfile;
	 // mfile.open("intpressure.txt", std::ios::app);
	 // mfile<<transition_depths.size()<<"  "<<ol_index<<std::endl;
	 // mfile.close();
        }
        else
          out.thermal_conductivities[i] = k_value;

        /*thermal expansivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
         For now only uses fperovskite periclase phase to calculate.*/
        if(thermal_alpha == 0.0)
          out.thermal_expansion_coefficients[i] = (a0+(a1*temperature)+a2*pow(temperature,-2))*exp(-a3*pressure*1e-9);
        else
          out.thermal_expansion_coefficients[i] = thermal_alpha;


        if(in.strain_rate.size())
        {
                      out.viscosities[i] = std::min(std::max(min_eta,calculate_viscosity(pressure,
                                                                       temperature,
                                                                       position,
                                                                       in.strain_rate[i])),max_eta);
        }


       double phase_dependence = 0.0;
       if (this->get_adiabatic_conditions().is_initialized())
        {
          for(unsigned int i=0; i<number_of_phase_transitions; ++i)
          {
            const double phaseFunction = phase_function (position,
                                                         temperature,
                                                         pressure,
                                                         i);

            phase_dependence += phaseFunction * density_jumps[i];
          }  
        }
        //density equation with pressure and temperature dependence, likely will change when adiabatic conditions are introduced.
          double density_temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            {
              // temperature dependence is 1 - alpha * (T - T(adiabatic))
              density_temperature_dependence -= (temperature - this->get_adiabatic_conditions().temperature(position))
                                                * out.thermal_expansion_coefficients[i];
            }
          else
              density_temperature_dependence -= temperature * out.thermal_expansion_coefficients[i];

          const double kappa = reference_compressibility;
          const double pressure_dependence = reference_rho * kappa * (pressure - this->get_surface_pressure());
          double rho = (reference_rho + pressure_dependence + phase_dependence)
                               * density_temperature_dependence;

          // in the end, all the influences are added up 
          out.densities[i] = rho;
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
            prm.declare_entry ("a0", "2.68e-5",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/K");
            prm.declare_entry ("a1", "2.77e-9",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/K^2");
            prm.declare_entry ("a2", "-1.21",
                               Patterns::Double(),
                               "coefficient for depth dependent thermal expansivity" "Units: K"); 
            prm.declare_entry ("a3", "3.76e-7",
                               Patterns::Double(0),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/GPa");
            prm.declare_entry ("c0", "3.48",
                               Patterns::List (Patterns::Double(0)),
                               "coefficient for depth dependent thermal conductivity" "Units: Wm^-1K^-1");
            prm.declare_entry ("c1", "5.17e-6",
                               Patterns::List (Patterns::Double(0)),
                               "coefficient for depth dependent thermal conductivity" "Units: Wm^-1K^-1GPa^-1");
            prm.declare_entry ("c2", "0.31",
                               Patterns::List (Patterns::Double(0)),
                               "coefficient for depth dependent thermal conductivity" "Units: None");

            //phase variables
            prm.declare_entry ("Auto transition temperature", "true",
                               Patterns::Bool (),
                               "If set to true the program will find the temperature at the transition "
                               "depth and use that for transition temperastures");
            prm.declare_entry ("Phase transition temperatures", "",
                               Patterns::List (Patterns::Double(0)),
                               "A list of temperatures where phase transitions occur. Higher or lower "
                               "temperatures lead to phase transition occurring in smaller or greater "
                               "depths than given in Phase transition depths, depending on the "
                               "Clapeyron slope given in Phase transition Clapeyron slopes. "
                               "List must have the same number of entries as Phase transition depths. "
                               "Units: $K$.");
            prm.declare_entry ("Phase transition Clapeyron slopes", "",
                               Patterns::List (Patterns::Double()),
                               "A list of Clapeyron slopes for each phase transition. A positive "
                               "Clapeyron slope indicates that the phase transition will occur in "
                               "a greater depth, if the temperature is higher than the one given in "
                               "Phase transition temperatures and in a smaller depth, if the "
                               "temperature is smaller than the one given in Phase transition temperatures. "
                               "For negative slopes the other way round. "
                               "List must have the same number of entries as Phase transition depths. "
                               "Units: $Pa/K$.");
            prm.declare_entry ("Phase transition density jumps", "",
                               Patterns::List (Patterns::Double(0)),
                               "A list of density jumps at each phase transition. A positive value means "
                               "that the density increases with depth. The corresponding entry in "
                               "Corresponding phase for density jump determines if the density jump occurs "
                               "in peridotite, eclogite or none of them."
                               "List must have the same number of entries as Phase transition depths. "
                               "Units: $kg/m^3$.");
            prm.declare_entry ("Phase transition widths", "",
                               Patterns::List (Patterns::Double(0)),
                               "A list of widths for each phase transition, in terms of depth. The phase functions "
                               "are scaled with these values, leading to a jump between phases "
                               "for a value of zero and a gradual transition for larger values. "
                               "List must have the same number of entries as Phase transition depths. "
                               "Units: $m$.");
            prm.declare_entry ("Phase transition depths", "",
                               Patterns::List (Patterns::Double(0)),
                               "A list of depths where phase transitions occur. Values must "
                               "monotonically increase. "
                               "Units: $m$.");


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
          a0                         = prm.get_double ("a0");
          a1                         = prm.get_double ("a1");
          a2                         = prm.get_double ("a2");
          a3                         = prm.get_double ("a3");
          c0                         = Utilities::string_to_double
                                       (Utilities::split_string_list(prm.get ("c0")));
          c1                         = Utilities::string_to_double
                                       (Utilities::split_string_list(prm.get ("c1")));
          c2                         = Utilities::string_to_double
                                       (Utilities::split_string_list(prm.get ("c2")));

          //phase transition parameters
          auto_temp                  = prm.get_bool ("Auto transition temperature");
          transition_depths          = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("Phase transition depths")));
          transition_widths          = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("Phase transition widths")));
          transition_temperatures    = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("Phase transition temperatures")));
          transition_slopes          = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("Phase transition Clapeyron slopes")));
          density_jumps              = Utilities::string_to_double
                                        (Utilities::split_string_list(prm.get ("Phase transition density jumps")));

              if (transition_widths.size() != transition_depths.size() ||
                  transition_temperatures.size() != transition_depths.size() ||
                  transition_slopes.size() != transition_depths.size() ||
                  density_jumps.size() != transition_depths.size())
                  AssertThrow(false, ExcMessage("Error: At least one list that gives input parameters for the phase "
                                              "transitions has the wrong size. Currently checking against transition depths. "
                                              "If phase transitions in terms of pressure inputs are desired, check to make sure "
                                              "'Define transition by depth instead of pressure = false'."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Declare dependencies on solution variables
      this->model_dependence.viscosity = NonlinearDependence::temperature | NonlinearDependence::pressure;
      this->model_dependence.density = NonlinearDependence::temperature | NonlinearDependence::pressure | NonlinearDependence::compositional_fields;
      this->model_dependence.compressibility = NonlinearDependence::none;
      this->model_dependence.specific_heat = NonlinearDependence::temperature | NonlinearDependence::pressure;
      this->model_dependence.thermal_conductivity = NonlinearDependence::temperature | NonlinearDependence::pressure;
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

