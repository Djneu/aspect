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
                          const SymmetricTensor<2,dim> &strain_rate,
                          const int phase) const
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
      double diffusion_energy_term = exp((diffusion_activation_energy[phase] + diffusion_activation_volume[phase] * adiabatic_pressure)
                         / (1 * constants::gas_constant * temperature));

      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double diffusion_adiabatic_energy_term
            = exp((diffusion_activation_energy[phase] + diffusion_activation_volume[phase] * adiabatic_pressure)
              / (1 * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double diffusion_temperature_dependence = diffusion_energy_term / diffusion_adiabatic_energy_term;
          if (diffusion_temperature_dependence > 1e6)
            diffusion_energy_term = diffusion_adiabatic_energy_term * 1e6;
          if (diffusion_temperature_dependence < 1.0 / 1e6)
            diffusion_energy_term = diffusion_adiabatic_energy_term / 1e6;
        }

      const double diffusion_strain_rate_dependence = (1.0 - 1) / 1;
      double vdiff = pow(diffusion_prefactor[phase],-1.0/1)
             * std::pow(second_strain_rate_invariant, diffusion_strain_rate_dependence)
             * pow(grain_size[phase], 3/1)
             * diffusion_energy_term;


      //energy term for calculation of viscosity from dislocation creep
      double dislocation_energy_term = exp((dislocation_activation_energy[phase] + dislocation_activation_volume[phase] * adiabatic_pressure)
                         / (3.5 * constants::gas_constant * temperature));
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double dislocation_adiabatic_energy_term
            = exp((dislocation_activation_energy[phase] + dislocation_activation_volume[phase] * adiabatic_pressure)
              / (3.5 * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double dislocation_temperature_dependence = dislocation_energy_term / dislocation_adiabatic_energy_term;
          if (dislocation_temperature_dependence > 1e6)
            dislocation_energy_term = dislocation_adiabatic_energy_term * 1e6;
          if (dislocation_temperature_dependence < 1.0 / 1e6)
            dislocation_energy_term = dislocation_adiabatic_energy_term / 1e6;
        }

      const double dislocation_strain_rate_dependence = (1.0 - 3.5) / 3.5;


      //if the strain rate inviariant is high enough calculate dislocation creep
      if(std::abs(second_strain_rate_invariant) > 1e-30)
      {
             vdis = pow(dislocation_prefactor[phase],-1.0/3.5)
             * std::pow(second_strain_rate_invariant,dislocation_strain_rate_dependence)
             * dislocation_energy_term;
      }

      //if dislocation creep was calculated use a harmonic average of the two, otherwise use diffusion creep.
      if(std::abs(second_strain_rate_invariant) > 1e-30)
        current_viscosity = vdis * vdiff / (vdis + vdiff);
      else
        current_viscosity = vdiff;

     return current_viscosity;
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    phase_function (const Point<dim> &position,
                    const double temperature,
                    const double pressure,
                    const int phase) const
    {

      double transition_pressure;
      double pressure_width;
      double width_temp;
      double phase_func;

      //if adiabatic profile has been created use that to calculate phase transitions based on pressure
      if(this->get_adiabatic_conditions().is_initialized())
        {
          //converting depth and width given to pressures
          const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
          const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + 0.5 * transition_widths[phase]);
          const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - 0.5 * transition_widths[phase]);
          transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
          pressure_width = (this->get_adiabatic_conditions().pressure(transition_plus_width)
                            - this->get_adiabatic_conditions().pressure(transition_minus_width));
          width_temp = transition_widths[phase];

          // then calculate the deviation from the transition point (both in temperature
          // and in pressure)
          double pressure_deviation = pressure - transition_pressure
                                      - transition_slopes[phase] * (temperature - transition_temperatures[phase]);

          // last, calculate the percentage of material that has undergone the transition
          // (also in dependence of the phase transition width - this is an input parameter)
          // use delta function for width = 0
          if (width_temp==0)
            (pressure_deviation > 0) ? phase_func = 1 : phase_func = 0;
          else
            phase_func = 0.5*(1.0 + std::tanh(pressure_deviation / pressure_width));
          return phase_func;
        }


          //for first time step use depth to calculate initial phases
      else
        {
          double depth = this->get_geometry_model().depth(position);
          double depth_deviation = (pressure > 0
                                    ?
                                    depth - transition_depths[phase]
                                    - transition_slopes[phase] * (depth / pressure) * (temperature - transition_temperatures[phase])
                                    :
                                    depth - transition_depths[phase]
                                    - transition_slopes[phase] / (this->get_gravity_model().gravity_vector(position).norm() * reference_rho)
                                    * (temperature - transition_temperatures[phase]));
          double phase_func;

          // use delta function for width = 0
          if (transition_widths[phase]==0)
            (depth_deviation > 0) ? phase_func = 1 : phase_func = 0;
          else
            phase_func = 0.5*(1.0 + std::tanh(depth_deviation / transition_widths[phase]));
          return phase_func;
        }
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    Pphase_function_derivative (const Point<dim> &,
                               const double temperature,
                               const double pressure,
                               unsigned int phase) const
    {
      double transition_pressure;
      double pressure_width;
      double width_temp;

      // we already should have the adiabatic conditions here
      AssertThrow (this->get_adiabatic_conditions().is_initialized(),
                   ExcMessage("need adiabatic conditions to incorporate phase transitions"));

      // first, get the pressure at which the phase transition occurs normally

      //converting depth and width given to pressures
      const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
      const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + 0.5 * transition_widths[phase]);
      const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - 0.5 * transition_widths[phase]);
      transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
      pressure_width = (this->get_adiabatic_conditions().pressure(transition_plus_width)
                        - this->get_adiabatic_conditions().pressure(transition_minus_width));
      width_temp = transition_widths[phase];

      // then calculate the deviation from the transition point (both in temperature
      // and in pressure)
      double pressure_deviation = pressure - transition_pressure
                                  - transition_slopes[phase] * (temperature - transition_temperatures[phase]);

      // last, calculate the analytical derivative of the phase function
      if (width_temp==0)
        return 0;
      else
        return 0.5 / pressure_width * (1.0 - std::tanh(pressure_deviation / pressure_width)
                                       * std::tanh(pressure_deviation / pressure_width));
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    Tphase_function_derivative (const Point<dim> &,
                               const double temperature,
                               const double pressure,
                               unsigned int phase) const
    {
      double transition_pressure;
      double pressure_width;
      double width_temp;

      // we already should have the adiabatic conditions here
      AssertThrow (this->get_adiabatic_conditions().is_initialized(),
                   ExcMessage("need adiabatic conditions to incorporate phase transitions"));

      // first, get the pressure at which the phase transition occurs normally

      //converting depth and width given to pressures
      const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
      const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + 0.5 * transition_widths[phase]);
      const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - 0.5 * transition_widths[phase]);
      transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
      pressure_width = (this->get_adiabatic_conditions().pressure(transition_plus_width)
                        - this->get_adiabatic_conditions().pressure(transition_minus_width));
      width_temp = transition_widths[phase];

      // then calculate the deviation from the transition point (both in temperature
      // and in pressure)
      double pressure_deviation = pressure - transition_pressure
                                  - transition_slopes[phase] * (temperature - transition_temperatures[phase]);
   
      // last, calculate the analytical derivative of the phase function
      if (width_temp==0)
        return 0;
      else
        return (( -0.5 * transition_slopes[phase] ) / pressure_width) * (1.0 - std::pow(std::tanh(pressure_deviation / pressure_width) , 2 ));
    }


    template <int dim>
    unsigned int
    PhaseTransitions<dim>::
    get_phase_index (const Point<dim> &position,
                     const double temperature,
                     const double pressure) const
    {
      unsigned int phase_index = 0;

      //Calculates for other phase transitions
      for(unsigned int j=1;j<transition_depths.size();++j)
        if(phase_function(position, temperature, pressure, j-1) > 0.5 )
          phase_index = j;

      //check to see if material has passed final phase transition
      if(transition_depths.size()>0)
        if(phase_function(position, temperature, pressure, transition_depths.size()-1) > 0.5)
          phase_index = transition_depths.size();


       //return which phase the material is in
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
        unsigned int ol_index = get_phase_index(position, temperature, pressure);
        double depth = this->get_geometry_model().depth(position);


        //constant properties
        out.compressibilities[i] = reference_compressibility;
        out.specific_heat[i] = reference_specific_heat;

        
        //Calculations for thermal conductivity and expansion coefficient
        double alpha;
        double conductivity;
        //thermal conductivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
        if(k_value == 0.0)
        {
          conductivity = (c0[ol_index]+(c1[ol_index]*pressure*1e-9))*pow((300/temperature),c2[ol_index]);
        }
        else
        {
          conductivity=k_value;
        }


        //thermal expansivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
       if(thermal_alpha == 0.0)
       {
        alpha = (a0[ol_index]+(a1[ol_index]*temperature)+a2[ol_index]*pow(temperature,-2))*exp(-a3[ol_index]*pressure*1e-9);

       }
       else
        {
          alpha = thermal_alpha;
        }
        out.thermal_expansion_coefficients[i] = alpha;
        out.thermal_conductivities[i] = conductivity;  


        //calculating viscosity
        double viscosity=0;
        if(in.strain_rate.size())
        {
                      viscosity = std::min(std::max(min_eta,calculate_viscosity(pressure,
                                                                       temperature,
                                                                       position,
                                                                       in.strain_rate[i],
                                                                       ol_index)),max_eta);
        }


   
        double density_phase_dependence = 0.0;
        double viscosity_phase_dependence = 1.0;
        for(unsigned int i=0; i<number_of_phase_transitions; ++i)
        {
          const double phaseFunction = phase_function (position,
                                                       temperature,
                                                       pressure,
                                                       i);

          density_phase_dependence += phaseFunction * density_jumps[i];
          viscosity_phase_dependence *= 1. + phaseFunction * (phase_prefactors[i+1]-1.);
        }  

        //density equation with pressure and temperature dependence, likely will change when adiabatic conditions are introduced.
         double density_phase_deviation = 0;
         if (this->get_adiabatic_conditions().is_initialized() && this->include_latent_heat())
           for (unsigned int ph=0; ph<number_of_phase_transitions; ++ph)
             {
                // calculate derivative of the phase function
                const double phase_derivative = Tphase_function_derivative(position,
                                                                       temperature,
                                                                       pressure,
                                                                       ph); 

               density_phase_deviation += phase_derivative*density_jumps[ph];
             }

          const double temperature_deviation = temperature - this->get_adiabatic_conditions().temperature(position);
          const double pressure_dev = pressure - this->get_adiabatic_conditions().pressure(position);
          double density_profile = reference_rho;
          if(this->get_adiabatic_conditions().is_initialized())
              density_profile = this->get_adiabatic_conditions().density(position);
        
          out.densities[i] = density_profile * (1 - alpha * temperature_deviation + reference_compressibility * pressure_dev)
                                  + density_phase_deviation * temperature_deviation ;

         //out.viscosities[i] = std::min(std::max(min_eta, viscosity * viscosity_phase_dependence), max_eta);
           out.viscosities[i] = viscosity;


        // Calculate entropy derivative
        {
          double entropy_gradient_pressure = 0.0;
          double entropy_gradient_temperature = 0.0;
          const double rho = out.densities[i];



          if (this->get_adiabatic_conditions().is_initialized() && this->include_latent_heat())
            for (unsigned int phase=0; phase<number_of_phase_transitions; ++phase)
              {
                // calculate derivative of the phase function
                const double PhaseFunctionDerivative = Pphase_function_derivative(position,
                                                                                 temperature,
                                                                                 pressure,
                                                                                 phase);

                // calculate the change of entropy across the phase transition
                double entropy_change = 0.0;
                entropy_change = transition_slopes[phase] * density_jumps[phase] / (rho * rho);


                // we need DeltaS * DX/Dpressure_deviation for the pressure derivative
                // and - DeltaS * DX/Dpressure_deviation * gamma for the temperature derivative
                entropy_gradient_pressure += PhaseFunctionDerivative * entropy_change;
                entropy_gradient_temperature -= PhaseFunctionDerivative * entropy_change * transition_slopes[phase];
              }
          out.entropy_derivative_pressure[i] = entropy_gradient_pressure;
          out.entropy_derivative_temperature[i] = entropy_gradient_temperature;
        }


        // Assign reaction terms
        for (unsigned int c=0; c<in.composition[i].size(); ++c)
          out.reaction_terms[i][c] = 0.0;

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

            // Grain size rheology parameters 
            prm.declare_entry ("Grain growth activation energy", "3.5e5",
                               Patterns::List (Patterns::Double(0)),
                               "The activation energy for grain growth $E_g$. "
                               "Units: $J/mol$.");
            prm.declare_entry ("Grain growth activation volume", "8e-6",
                               Patterns::List (Patterns::Double(0)),
                               "The activation volume for grain growth $E_g$. "
                               "Units: $m^3/mol$.");
            prm.declare_entry ("Grain growth exponent", "3",
                               Patterns::List (Patterns::Double(0)),
                               "Exponent of the grain growth law $p_g$. This is an experimentally determined "
                               "grain growth constant. "
                               "Units: none.");
            prm.declare_entry ("Grain growth rate constant", "1.5e-5",
                               Patterns::List (Patterns::Double(0)),
                               "Prefactor of the Ostwald ripening grain growth law $G_0$. "
                               "This is dependent on water content, which is assumed to be "
                               "50 H/10^6 Si for the default value. "
                               "Units: $m^{p_g}/s$.");
            prm.declare_entry ("Reciprocal required strain", "10",
                               Patterns::List (Patterns::Double(0)),
                               "This parameters $\\lambda$ gives an estimate of the strain necessary "
                               "to achieve a new grain size. ");
            prm.declare_entry ("Recrystallized grain size", "",
                               Patterns::List (Patterns::Double(0)),
                               "The grain size $d_{ph}$ to that a phase will be reduced to when crossing a phase transition. "
                               "When set to zero, grain size will not be reduced. "
                               "Units: m.");
            prm.declare_entry ("Use paleowattmeter", "true",
                               Patterns::Bool (),
                               "A flag indicating whether the computation should be use the "
                               "paleowattmeter approach of Austin and Evans (2007) for grain size reduction "
                               "in the dislocation creep regime (if true) or the paleopiezometer aprroach "
                               "from Hall and Parmetier (2003) (if false).");
            prm.declare_entry ("Average specific grain boundary energy", "1.0",
                               Patterns::List (Patterns::Double(0)),
                               "The average specific grain boundary energy $\\gamma$. "
                               "Units: J/m^2.");
            prm.declare_entry ("Work fraction for boundary area change", "0.1",
                               Patterns::List (Patterns::Double(0)),
                               "The fraction $\\chi$ of work done by dislocation creep to change the grain boundary area. "
                               "Units: J/m^2.");
            prm.declare_entry ("Geometric constant", "3",
                               Patterns::List (Patterns::Double(0)),
                               "Geometric constant $c$ used in the paleowattmeter grain size reduction law. "
                               "Units: none.");
            prm.declare_entry ("Grain size", "0.003",Patterns::List (Patterns::Double(0)),
                               "Stabilizes strain dependent viscosity. Units: m");
          prm.declare_entry ("Minimum grain size", "1e-5",
                             Patterns::Double (0),
                             "The minimum grain size that is used for the material model. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allows for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: Pa s.");
          prm.declare_entry ("Lower mantle grain size scaling", "1.0",
                             Patterns::Double (0),
                             "A scaling factor for the grain size in the lower mantle. In models where the "
                             "high grain size contrast between the upper and lower mantle causes numerical "
                             "problems, the grain size in the lower mantle can be scaled to a larger value, "
                             "simultaneously scaling the viscosity prefactors and grain growth parameters "
                             "to keep the same physical behavior. Differences to the original formulation "
                             "only occur when material with a smaller grain size than the recrystallization "
                             "grain size cross the upper-lower mantle boundary. "
                             "The real grain size can be obtained by dividing the model grain size by this value. "
                             "Units: none.");
          prm.declare_entry ("Advect logarithm of grain size", "false",
                             Patterns::Bool (),
                             "Whether to advect the logarithm of the grain size or the "
                             "grain size. The equation and the physics are the same, "
                             "but for problems with high grain size gradients it might "
                             "be preferable to advect the logarithm. ");


            // Discloation viscosity parameters
            prm.declare_entry ("Dislocation viscosity iteration threshold", "1e-3",
                               Patterns::Double(0),
                               "We need to perform an iteration inside the computation "
                               "of the dislocation viscosity, because it depends on the "
                               "dislocation strain rate, which depends on the dislocation "
                               "viscosity itself. This number determines the termination "
                               "accuracy, i.e. if the dislocation viscosity changes by less "
                               "than this factor we terminate the iteration.");
            prm.declare_entry ("Dislocation viscosity iteration number", "10",
                               Patterns::Integer(0),
                               "We need to perform an iteration inside the computation "
                               "of the dislocation viscosity, because it depends on the "
                               "dislocation strain rate, which depends on the dislocation "
                               "viscosity itself. This number determines the maximum "
                               "number of iterations that are performed. ");
            prm.declare_entry ("Dislocation prefactor", "1.57e-017",
                               Patterns::List (Patterns::Double(0)),
                               "Discloation viscosity prefactor" "Units: None");
            prm.declare_entry ("Dislocation activation energy", "530000",
                               Patterns::List (Patterns::Double(0)),
                               "Activation energy for dislocation viscosity." "Units: J/mol");
            prm.declare_entry ("Dislocation activation volume", "1.4e-005",
                               Patterns::List (Patterns::Double(0)),
                               "Activation volume for dislocation viscosity." "Units: m^3/mol");
           prm.declare_entry ("Dislocation creep exponent", "3.5",
                             Patterns::List (Patterns::Double(0)),
                             "Power-law exponent $n_{dis}$ for dislocation creep. "
                             "Units: none.");

            //Diffusion viscosity parameters
            prm.declare_entry ("Diffusion prefactor", "1.25e-015",
                               Patterns::List (Patterns::Double(0)),
                               "Diffusion viscosity prefactor with first entry being diffusion, second dislocation" "Units: None");
            prm.declare_entry ("Diffusion activation energy", "375000",
                               Patterns::List (Patterns::Double(0)),
                               "Activation energy for viscosity equation." "Units: J/mol");
            prm.declare_entry ("Diffusion activation volume", "6e-6",
                               Patterns::List (Patterns::Double(0)),
                               "Diffusion activation volume for viscosity equation." "Units: m^3/mol");
          prm.declare_entry ("Diffusion creep grain size exponent", "3",
                             Patterns::List (Patterns::Double(0)),
                             "Diffusion creep grain size exponent $p_{diff}$ that determines the "
                             "dependence of vescosity on grain size. "
                             "Units: none.");

            // Additional viscosity parameters
          prm.declare_entry ("Maximum temperature dependence of viscosity", "100",
                             Patterns::Double (0),
                             "The factor by which viscosity at adiabatic temperature and ambient temperature "
                             "are allowed to differ (a value of x means that the viscosity can be x times higher "
                             "or x times lower compared to the value at adiabatic temperature. This parameter "
                             "is introduced to limit local viscosity contrasts, but still allow for a widely "
                             "varying viscosity over the whole mantle range. "
                             "Units: none.");
            prm.declare_entry ("Maximum viscosity","1.0e24",Patterns::Double(0),
                               "Reference strain rate for first time step. Units: Kg/m/s");
            prm.declare_entry ("Minimum viscosity", "1.0e18", Patterns::Double(0),
                               "Stabilizes strain dependent viscosity. Units: Kg/m/s");

            //thermal expansivity and conductivity parameters
            prm.declare_entry ("a0", "2.68e-5",
                               Patterns::List (Patterns::Double(0)),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/K");
            prm.declare_entry ("a1", "2.77e-9",
                               Patterns::List (Patterns::Double(0)),
                               "coefficient for depth dependent thermal expansivity" "Units: 1/K^2");
            prm.declare_entry ("a2", "-1.21",
                               Patterns::List (Patterns::Double()),
                               "coefficient for depth dependent thermal expansivity" "Units: K"); 
            prm.declare_entry ("a3", "3.76e-7",
                               Patterns::List (Patterns::Double(0)),
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
            prm.declare_entry ("Phase transition temperature jumps", "",
                               Patterns::List (Patterns::Double()),
                               "A list of temperature jumps at each phase transition. A positive value means "
                               "that the temperature increases with depth. The corresponding entry in "
                               "Corresponding phase for temperature jump determines if the temperature jump occurs "
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
            prm.declare_entry ("Viscosity prefactors", "",
                               Patterns::List (Patterns::Double(0)),
                               "A list of prefactors for the viscosity for each phase. The reference "
                               "viscosity will be multiplied by this factor to get the corresponding "
                               "viscosity for each phase. "
                               "List must have one more entry than Phase transition depths. "
                               "Units: non-dimensional.");


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

          // grain evolution parameters
          grain_growth_activation_energy        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation energy")));
          grain_growth_activation_volume        = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth activation volume")));
          grain_growth_rate_constant            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth rate constant")));
          grain_growth_exponent                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Grain growth exponent")));
          reciprocal_required_strain            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Reciprocal required strain")));
          min_grain_size                        = prm.get_double ("Minimum grain size");
          pv_grain_size_scaling                 = prm.get_double ("Lower mantle grain size scaling");
          use_paleowattmeter                    = prm.get_bool ("Use paleowattmeter");
          grain_boundary_energy                 = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Average specific grain boundary energy")));
          boundary_area_change_work_fraction    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Work fraction for boundary area change")));
          geometric_constant                    = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Geometric constant")));
          advect_log_gransize                   = prm.get_bool ("Advect logarithm of grain size");

          //viscosity parameters
          dislocation_viscosity_iteration_threshold = prm.get_double("Dislocation viscosity iteration threshold");
          dislocation_viscosity_iteration_number = prm.get_integer("Dislocation viscosity iteration number");
          max_temperature_dependence_of_eta     = prm.get_double ("Maximum temperature dependence of viscosity");
          max_eta                              = prm.get_double ("Maximum viscosity");
          min_eta                              = prm.get_double("Minimum viscosity");
          grain_size                           = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Grain size")));
          diffusion_prefactor                  = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Diffusion prefactor")));
          diffusion_activation_energy          = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Diffusion activation energy")));
          diffusion_activation_volume          = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Diffusion activation volume")));
          diffusion_creep_grain_size_exponent   = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Diffusion creep grain size exponent")));
          dislocation_prefactor                = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Dislocation prefactor")));
          dislocation_activation_energy        = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Dislocation activation energy")));
          dislocation_activation_volume        = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Dislocation activation volume")));
          dislocation_creep_exponent            = Utilities::string_to_double
                                                  (Utilities::split_string_list(prm.get ("Dislocation creep exponent")));

          if(diffusion_prefactor.size() != diffusion_activation_energy.size() ||
             diffusion_prefactor.size() != diffusion_activation_volume.size() || 
             diffusion_prefactor.size() != dislocation_prefactor.size() || 
             diffusion_prefactor.size() != dislocation_activation_energy.size() || 
             diffusion_prefactor.size() != dislocation_activation_volume.size())
               AssertThrow(false, ExcMessage("Error: Too many inputs for activation volume, activation energy, or viscosity prefactor."));

          //thermal exansivity and conductivity parameters
          a0                         = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("a0")));
          a1                         = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("a1")));
          a2                         = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("a2")));
          a3                         = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("a3")));
          c0                         = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("c0")));
          c1                         = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("c1")));
          c2                         = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("c2")));

              if (a0.size() != a1.size() ||
                  a1.size() != a2.size() ||
                  a2.size() != a3.size() ||
                  a0.size() != c0.size() ||  
                  c0.size() != c1.size() ||
                  c1.size() != c2.size())
                  AssertThrow(false, ExcMessage("Error, one of your a or c values is not right"));

          //phase transition parameters
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
          temperature_jumps              = Utilities::string_to_double
                                        (Utilities::split_string_list(prm.get ("Phase transition temperature jumps")));
          phase_prefactors            = Utilities::string_to_double
                                          (Utilities::split_string_list(prm.get ("Viscosity prefactors")));

              if (transition_widths.size() != transition_depths.size() ||
                  transition_temperatures.size() != transition_depths.size() ||
                  transition_slopes.size() != transition_depths.size() ||
                  density_jumps.size() != transition_depths.size() ||
                  phase_prefactors.size() != density_jumps.size()+1)
                  AssertThrow(false, ExcMessage("Error: At least one list that gives input parameters for the phase "
                                              "transitions has the wrong size. Currently checking against transition depths. "
                                              "If phase transitions in terms of pressure inputs are desired, check to make sure "
                                              "'Define transition by depth instead of pressure = false'."));

          // as the phase viscosity prefactors are all applied multiplicatively on top of each other,
          // we have to scale them here so that they are relative factors in comparison to the product
          // of the prefactors of all phase above the current one
          for (unsigned int phase=1; phase<phase_prefactors.size(); ++phase)
            {
              phase_prefactors[phase] /= phase_prefactors[phase-1];
            }
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

