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
#include <deal.II/base/signaling_nan.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/newton.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace MaterialModel
  {

  namespace
   {
     std::vector<std::string> make_phase_additional_outputs_names()
     {
       std::vector<std::string> names;
       names.emplace_back("phase");
       names.emplace_back("diffusion");
       names.emplace_back("dislocation");
       names.emplace_back("viscosity_ratio");
       return names;
     }
   }

   template <int dim>
   PhaseAdditionalOutputs<dim>::PhaseAdditionalOutputs (const unsigned int n_points)
     :
     NamedAdditionalMaterialOutputs<dim>(make_phase_additional_outputs_names()),
     phase(n_points, numbers::signaling_nan<double>()),
     diffusion(n_points, numbers::signaling_nan<double>()),
     dislocation(n_points, numbers::signaling_nan<double>()),
     viscosity_ratio(n_points, numbers::signaling_nan<double>())
   {}

   template <int dim>
   std::vector<double>
   PhaseAdditionalOutputs<dim>::get_nth_output(const unsigned int idx) const
   {
     AssertIndexRange (idx, 4);
     switch (idx)
       {
         case 0:
    	   return phase;

         case 1:
           return diffusion;

         case 2:
           return dislocation;

         case 3:
           return viscosity_ratio;

         default:
           AssertThrow(false, ExcInternalError());
       }
     // We will never get here, so just return something
     return phase;
   }



    template <int dim>
    void
    PhaseTransitions<dim>::
    convert_log_grain_size (std::vector<double> &composition) const
    {
        // get grain size and limit it to a global minimum
        const unsigned int grain_size_index = this->introspection().compositional_index_for_name("grain_size");
        double grain_size = composition[grain_size_index];
        grain_size = std::max(std::exp(-grain_size),min_grain_size);

        composition[grain_size_index] = grain_size;
    }


    template <int dim>
    double
    PhaseTransitions<dim>::
    grain_size_growth_rate (const double                  temperature,
                            const double                  pressure,
                            const std::vector<double>    &compositional_fields,
                            const SymmetricTensor<2,dim> &strain_rate,
                            const Tensor<1,dim>          &/*velocity*/,
                            const Point<dim>             &position,
                            const unsigned int            field_index,
                            const int                     crossed_transition) const
    {
        // we want to iterate over the grain size evolution here, as we solve in fact an ordinary differential equation
        // and it is not correct to use the starting grain size (and introduces instabilities)
        const double original_grain_size = compositional_fields[field_index];
        if ((original_grain_size != original_grain_size) || this->get_timestep() == 0.0
            || original_grain_size < std::numeric_limits<double>::min())
          return 0.0;

        // set up the parameters for the sub-timestepping of grain size evolution
        std::vector<double> current_composition = compositional_fields;
        double grain_size = original_grain_size;
        double grain_size_change = 0.0;
        const double timestep = this->get_timestep();

        // use a sub timestep of 500 yrs, currently fixed timestep
        double grain_growth_timestep = 500 * 3600 * 24 * 365.25;
        double time = 0;

        // find out in which phase we are
        const unsigned int phase_index = get_phase_index(position, temperature, pressure);

        // we keep the dislocation viscosity of the last iteration as guess
        // for the next one
        double current_dislocation_viscosity = 0.0;

        do
          {
            time += grain_growth_timestep;

            if (timestep - time < 0)
              {
                grain_growth_timestep = timestep - (time - grain_growth_timestep);
                time = timestep;
              }

            // grain size growth due to Ostwald ripening
            const double m = grain_growth_exponent[phase_index];
            const double grain_size_growth_rate = grain_growth_rate_constant[phase_index] / (m * pow(grain_size,m-1))
                                                  * exp(- (grain_growth_activation_energy[phase_index] + pressure * grain_growth_activation_volume[phase_index])
                                                        / (constants::gas_constant * temperature));
            const double grain_size_growth = grain_size_growth_rate * grain_growth_timestep;

            // grain size reduction in dislocation creep regime
            const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
            const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

            const double current_diffusion_viscosity   = diffusion_viscosity(temperature, pressure, current_composition, strain_rate, position);
            current_dislocation_viscosity              = dislocation_viscosity(temperature, pressure, current_composition, strain_rate, position, current_dislocation_viscosity);

            double current_viscosity;
            if (std::abs(second_strain_rate_invariant) > 1e-30)
              current_viscosity = current_dislocation_viscosity * current_diffusion_viscosity / (current_dislocation_viscosity + current_diffusion_viscosity);
            else
              current_viscosity = current_diffusion_viscosity;

            const double dislocation_strain_rate = second_strain_rate_invariant
                                                   * current_viscosity / current_dislocation_viscosity;

            double grain_size_reduction = 0.0;

            if (use_paleowattmeter)
              {
                // paleowattmeter: Austin and Evans (2007): Paleowattmeters: A scaling relation for dynamically recrystallized grain size. Geology 35, 343-346
                const double stress = 2.0 * second_strain_rate_invariant * current_viscosity;
                const double grain_size_reduction_rate = 2.0 * stress * boundary_area_change_work_fraction[phase_index] * dislocation_strain_rate * pow(grain_size,2)
                                                         / (geometric_constant[phase_index] * grain_boundary_energy[phase_index]);
                grain_size_reduction = grain_size_reduction_rate * grain_growth_timestep;
              }
            else
              {
                // paleopiezometer: Hall and Parmentier (2003): Influence of grain size evolution on convective instability. Geochem. Geophys. Geosyst., 4(3).
                grain_size_reduction = reciprocal_required_strain[phase_index] * dislocation_strain_rate * grain_size * grain_growth_timestep;
              }

            grain_size_change = grain_size_growth - grain_size_reduction;

            // If the change in grain size is very large or small decrease timestep and try
            // again, or increase timestep and move on.
            if ((grain_size_change / grain_size < 0.001 && grain_size_growth / grain_size < 0.1
                 && grain_size_reduction / grain_size < 0.1) || grain_size == 0.0)
              grain_growth_timestep *= 2;
            else if (grain_size_change / grain_size > 0.1 || grain_size_growth / grain_size > 0.5
                     || grain_size_reduction / grain_size > 0.5)
              {
                grain_size_change = 0.0;
                time -= grain_growth_timestep;

                grain_growth_timestep /= 2.0;
              }

            //stop

            grain_size += grain_size_change;
            current_composition[field_index] = grain_size;

            Assert(grain_size > 0,
                   ExcMessage("The grain size became smaller than zero. This is not valid, "
                              "and likely an effect of a too large sub-timestep, or unrealistic "
                              "input parameters."));
          }
        while (time < timestep);

        // reduce grain size to recrystallized_grain_size when crossing phase transitions
        // if the distance in radial direction a grain moved compared to the last time step
        // is crossing a phase transition, reduce grain size

        // TODO: recrystallize first, and then do grain size growth/reduction for grains that crossed the transition
        // in dependence of the distance they have moved
        double phase_grain_size_reduction = 0.0;
        if (this->introspection().name_for_compositional_index(field_index) == "grain_size"
            &&
            this->get_timestep_number() > 0)
          {
            // check if material has crossed any phase transition, if yes, reset grain size
            if (crossed_transition != -1)
              if (recrystallized_grain_size[crossed_transition] > 0.0)
                phase_grain_size_reduction = grain_size - recrystallized_grain_size[crossed_transition];
          }

        grain_size = std::max(grain_size, min_grain_size);

        return grain_size - original_grain_size - phase_grain_size_reduction;
    }


    template <int dim>
    double
    PhaseTransitions<dim>::
    diffusion_viscosity (const double                  temperature,
                         const double                  pressure,
                         const std::vector<double>    &composition,
                         const SymmetricTensor<2,dim> &strain_rate,
                         const Point<dim>             &position) const
    {
      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      // find out in which phase we are
      const unsigned int ol_index = get_phase_index(position, temperature, pressure);

      // TODO: make this more general, for more phases we have to average grain size somehow
      // TODO: default when field is not given & warning
      // limit the grain size to a global minimum
      const std::string field_name = "grain_size";
      const double grain_size = this->introspection().compositional_name_exists(field_name)
                                ?
                                composition[this->introspection().compositional_index_for_name(field_name)]
                                :
                                constant_grain_size[ol_index];

      // TODO: we use the prefactors from Behn et al., 2009 as default values, but their laws use the strain rate
      // and we use the second invariant --> check if the prefactors should be changed
      double energy_term = exp((diffusion_activation_energy[ol_index] + diffusion_activation_volume[ol_index] * pressure)
                         / (1.0 * constants::gas_constant * temperature));

      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = exp((diffusion_activation_energy[ol_index] + diffusion_activation_volume[ol_index] * pressure)
              / (1.0 * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      //const double strain_rate_dependence = (1.0 - diffusion_creep_exponent[ol_index]) / diffusion_creep_exponent[ol_index];
      const double strain_rate_dependence = 0.;

      return pow(diffusion_prefactor[ol_index],-1.0/1.0)
             * std::pow(second_strain_rate_invariant,strain_rate_dependence)
             * pow(grain_size, diffusion_creep_grain_size_exponent[ol_index]/1.0)
             * energy_term;
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    dislocation_viscosity (const double      temperature,
                           const double      pressure,
                           const std::vector<double> &composition,
                           const SymmetricTensor<2,dim> &strain_rate,
                           const Point<dim> &position,
                           const double viscosity_guess) const
    {
      const double diff_viscosity = diffusion_viscosity(temperature,pressure,composition,strain_rate,position) ;

      // Start the iteration with the full strain rate
      double dis_viscosity;
      if (viscosity_guess == 0)
        dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,pressure,std::vector<double>(),strain_rate,position);
      else
        dis_viscosity = viscosity_guess;

      double dis_viscosity_old = 0;
      unsigned int i = 0;
      while (std::abs((dis_viscosity-dis_viscosity_old) / dis_viscosity) > dislocation_viscosity_iteration_threshold && i < dislocation_viscosity_iteration_number)
          {
          SymmetricTensor<2,dim> dislocation_strain_rate = diff_viscosity
                                             / (diff_viscosity + dis_viscosity) * strain_rate;
          dis_viscosity_old = dis_viscosity;
          dis_viscosity = dislocation_viscosity_fixed_strain_rate(temperature,
                                                                  pressure,
                                                                  std::vector<double>(),
                                                                  dislocation_strain_rate,
                                                                  position);
          i++;
          }
      return dis_viscosity;
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    viscosity (const double temperature,
               const double pressure,
               const std::vector<double> &composition,
               const SymmetricTensor<2,dim> &strain_rate,
               const Point<dim> &position) const
    {
      //TODO: add assert
      /*if (this->get_timestep_number() > 0)
        Assert (grain_size >= 1.e-6, ExcMessage ("Error: The grain size should not be smaller than 1e-6 m."));*/

      const SymmetricTensor<2,dim> shear_strain_rate = strain_rate - 1./dim * trace(strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      const double diff_viscosity = diffusion_viscosity(temperature, pressure, composition, strain_rate, position);

      double effective_viscosity;
      if(std::abs(second_strain_rate_invariant) > 1e-30 && use_dislocation == true)
        {
          const double disl_viscosity = dislocation_viscosity(temperature, pressure, composition, strain_rate, position);
          effective_viscosity = disl_viscosity * diff_viscosity / (disl_viscosity + diff_viscosity);
        }
      else
        effective_viscosity = diff_viscosity;

      return effective_viscosity;
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    dislocation_viscosity_fixed_strain_rate (const double      temperature,
                                             const double      pressure,
                                             const std::vector<double> &,
                                             const SymmetricTensor<2,dim> &dislocation_strain_rate,
                                             const Point<dim> &position) const
    {
      const SymmetricTensor<2,dim> shear_strain_rate = dislocation_strain_rate - 1./dim * trace(dislocation_strain_rate) * unit_symmetric_tensor<dim>();
      const double second_strain_rate_invariant = std::sqrt(std::abs(second_invariant(shear_strain_rate)));

      // Currently this will never be called without adiabatic_conditions initialized, but just in case
      const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                        ?
                                        this->get_adiabatic_conditions().pressure(position)
                                        :
                                        pressure;

      // find out in which phase we are
      const unsigned int ol_index = get_phase_index(position, temperature, adiabatic_pressure);

      double energy_term = exp((dislocation_activation_energy[ol_index] + dislocation_activation_volume[ol_index] * adiabatic_pressure)
                         / (dislocation_creep_exponent[ol_index] * constants::gas_constant * temperature));
      if (this->get_adiabatic_conditions().is_initialized())
        {
          const double adiabatic_energy_term
            = exp((dislocation_activation_energy[ol_index] + dislocation_activation_volume[ol_index] * adiabatic_pressure)
              / (dislocation_creep_exponent[ol_index] * constants::gas_constant * this->get_adiabatic_conditions().temperature(position)));

          const double temperature_dependence = energy_term / adiabatic_energy_term;
          if (temperature_dependence > max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term * max_temperature_dependence_of_eta;
          if (temperature_dependence < 1.0 / max_temperature_dependence_of_eta)
            energy_term = adiabatic_energy_term / max_temperature_dependence_of_eta;
        }

      const double strain_rate_dependence = (1.0 - dislocation_creep_exponent[ol_index]) / dislocation_creep_exponent[ol_index];

      return pow(dislocation_prefactor[ol_index],-1.0/dislocation_creep_exponent[ol_index])
             * std::pow(second_strain_rate_invariant,strain_rate_dependence)
             * energy_term;
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
    unsigned int
    PhaseTransitions<dim>::
    get_phase_index (const Point<dim> &position,
                     const double temperature,
                     const double pressure) const
    {
      unsigned int phase_index = 0;

      //Calculates for other phase transitions
      for(unsigned int j=1;j<transition_depths.size();++j)
        if(phase_function(position, temperature, pressure, j-1) >= 0.5 )
          phase_index = j;

      //check to see if material has passed final phase transition
      if(transition_depths.size()>0)
        if(phase_function(position, temperature, pressure, transition_depths.size()-1) >= 0.5)
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
        const SymmetricTensor<2,dim> strain_rate = in.strain_rate[i];

        //use adiabatic pressure for phase index
        const double adiabatic_pressure = this->get_adiabatic_conditions().is_initialized()
                                        ?
                                        this->get_adiabatic_conditions().pressure(position)
                                        :
                                        pressure;

        unsigned int ol_index = get_phase_index(position, temperature, adiabatic_pressure);

        // Reset entropy derivatives
        out.entropy_derivative_pressure[i] = 0;
        out.entropy_derivative_temperature[i] = 0;

        // Reset all reaction terms
        for (unsigned int c=0; c < in.composition[i].size(); ++c )
          {
            out.reaction_terms[i][c]            = 0;
          }

    	// convert the grain size from log to normal
    	std::vector<double> composition (in.composition[i]);
       /* if (advect_log_grainsize)
          convert_log_grain_size(composition);
        else
          {
            const unsigned int grain_size_index = this->introspection().compositional_index_for_name("grain_size");
            composition[grain_size_index] = std::max(min_grain_size,composition[grain_size_index]);
          }*/

        // set up an integer that tells us which phase transition has been crossed inside of the cell
        int crossed_transition(-1);

        if (this->get_adiabatic_conditions().is_initialized())
          for (unsigned int phase=0;phase<transition_depths.size();++phase)
            {
              // first, get the pressure at which the phase transition occurs normally
              const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[phase]);
              const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[phase] + transition_widths[phase]);
              const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[phase] - transition_widths[phase]);
              const double transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
              const double pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                                   - this->get_adiabatic_conditions().pressure(transition_minus_width));


              // then calculate the deviation from the transition point (both in temperature
              // and in pressure)
              double pressure_deviation = in.pressure[i] - transition_pressure
                                          - transition_slopes[phase] * (in.temperature[i] - transition_temperatures[phase]);

              if ((std::abs(pressure_deviation) < pressure_width)
                &&
                ((in.velocity[i] * this->get_gravity_model().gravity_vector(in.position[i])) * pressure_deviation > 0))
                crossed_transition = phase;
            }
        else
          for (unsigned int j=0; j<in.position.size(); ++j)
            for (unsigned int k=0;k<transition_depths.size();++k)
              if((phase_function(in.position[i], in.temperature[i], in.pressure[i], k)
                  != phase_function(in.position[j], in.temperature[j], in.pressure[j], k))
                  &&
                  ((in.velocity[i] * this->get_gravity_model().gravity_vector(in.position[i]))
                  * ((in.position[i] - in.position[j]) * this->get_gravity_model().gravity_vector(in.position[i])) > 0))
                crossed_transition = k;


        if (in.strain_rate.size() > 0)
        {
          out.viscosities[i] = std::min(max_eta, std::max(min_eta, viscosity(temperature,
                                  adiabatic_pressure,
                                  composition,
                                  strain_rate,
                                  position)));
        }



        
        //Calculations for thermal conductivity and expansion coefficient
        double alpha;
        double conductivity;
        //thermal conductivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
        if(k_value == 0.0)
          conductivity = (c0[ol_index]+(c1[ol_index]*adiabatic_pressure*1e-9))*pow((300/temperature),c2[ol_index]);
        else
          conductivity=k_value;


        //thermal expansivity equation (Tosi, et. al., 2013) that varies with temperature, depth, and phase.
       if(thermal_alpha == 0.0)
          alpha = (a0[ol_index]+(a1[ol_index]*temperature)+a2[ol_index]*pow(temperature,-2))*exp(-a3[ol_index]*adiabatic_pressure*1e-9);
       else
          alpha = thermal_alpha;


        double density_phase_dependence = 0.0;
        double viscosity_phase_dependence = 1.0;
        for(unsigned int i=0; i<number_of_phase_transitions; ++i)
        {
          const double phaseFunction = phase_function (position,
                                                       temperature,
                                                       adiabatic_pressure,
                                                       i);

          density_phase_dependence += phaseFunction * density_jumps[i];
          
          viscosity_phase_dependence *= 1. + phaseFunction * (phase_prefactors[i+1]-1.);
        }

          const double temperature_deviation = temperature - this->get_adiabatic_conditions().temperature(position);
          const double pressure_dev = pressure - this->get_adiabatic_conditions().pressure(position);
          double density_profile = reference_rho;
          if(this->get_adiabatic_conditions().is_initialized())
              density_profile = this->get_adiabatic_conditions().density(position);
        
          out.densities[i] = (density_profile) * (1 - alpha * temperature_deviation + reference_compressibility * pressure_dev) + density_phase_dependence;
                                  //tdensity_phase_deviation * temperature_deviation

          out.thermal_expansion_coefficients[i] = alpha;
          out.thermal_conductivities[i] = conductivity;
  
          //constant properties
          out.compressibilities[i] = reference_compressibility;
          out.specific_heat[i] = reference_specific_heat;


        // Calculate entropy derivative for phase_changes
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
                                                                                 adiabatic_pressure,
                                                                                 phase);

                // calculate the change of entropy across the phase transition
                double entropy_change = 0.0;
                entropy_change = transition_slopes[phase] * density_jumps[phase] / (rho * rho);


                // we need DeltaS * DX/Dpressure_deviation for the pressure derivative
                // and - DeltaS * DX/Dpressure_deviation * gamma for the temperature derivative
                entropy_gradient_pressure += PhaseFunctionDerivative * entropy_change;
                entropy_gradient_temperature -= PhaseFunctionDerivative * entropy_change * transition_slopes[phase];
              }

          // Entropy derivatives for phase changes
          out.entropy_derivative_pressure[i] = entropy_gradient_pressure;
          out.entropy_derivative_temperature[i] = entropy_gradient_temperature;

          // Melting 
          out.entropy_derivative_pressure[i] += entropy_derivative_melt (temperature, pressure, composition,
                                                                   position, NonlinearDependence::pressure) ; // for pressure dependence

          out.entropy_derivative_temperature[i] += entropy_derivative_melt (temperature, pressure, composition,
                                                                      position, NonlinearDependence::temperature) ; // for temperature dependence


        }

        // TODO: make this more general for not just olivine grains
        if (in.strain_rate.size() > 0)
          for (unsigned int c=0;c<composition.size();++c)
            {
             /* if (this->introspection().name_for_compositional_index(c) == "grain_size")
              {
                out.reaction_terms[i][c] = grain_size_growth_rate(in.temperature[i], in.pressure[i], composition,
                    in.strain_rate[i], in.velocity[i], in.position[i], c, crossed_transition);
                if(advect_log_grainsize)
                  out.reaction_terms[i][c] = - out.reaction_terms[i][c] / composition[c];
              }*/
            /*  if (this->introspection().name_for_compositional_index(c) == "grain_size")
              {
            	  //stop here
                  out.reaction_terms[i][c] = grain_size_growth_rate(in.temperature[i], in.pressure[i], composition,
                                                               in.strain_rate[i], in.velocity[i], in.position[i], c, crossed_transition);
                  if (advect_log_grainsize)
                    out.reaction_terms[i][c] = - out.reaction_terms[i][c] / composition[c];
              }*/
              if (this->introspection().name_for_compositional_index(c) == "peridotite_melt_fraction")
                {
                  out.reaction_terms[i][c] = peridotite_melt_fraction(in.temperature[i], in.pressure[i], composition, in.position[i]) - in.composition[i][c];
                }
              else if (this->introspection().name_for_compositional_index(c) == "test")
                {
                  out.reaction_terms[i][c] = 1;
                }
              else
                out.reaction_terms[i][c] = 0.0;

            }

          if (PhaseAdditionalOutputs<dim> *phase_out = out.template get_additional_output<PhaseAdditionalOutputs<dim> >())
            {
              phase_out->phase[i] = ol_index;
              phase_out->diffusion[i] = diffusion_viscosity(in.temperature[i],adiabatic_pressure,composition,in.strain_rate[i],in.position[i]);
              phase_out->dislocation[i] = dislocation_viscosity(in.temperature[i], adiabatic_pressure, composition, in.strain_rate[i], in.position[i]);
              phase_out->viscosity_ratio[i] = phase_out->diffusion[i]/phase_out->dislocation[i];
            }


      }
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    entropy_derivative_melt (const double temperature,
                             const double pressure,
                             const std::vector<double> &compositional_fields,
                             const Point<dim> &position,
                             const NonlinearDependence::Dependence dependence) const
    {
      double entropy_gradient = 0.0;

      // calculate latent heat of melting
      // we need the change of melt fraction in dependence of pressure and temperature

      // for peridotite after Katz, 2003
      const double T_solidus        = A1 + 273.15
                                      + A2 * pressure
                                      + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus       = C1 + 273.15
                                      + C2 * pressure
                                      + C3 * pressure * pressure;

      const double dT_solidus_dp        = A2 + 2 * A3 * pressure;
      const double dT_lherz_liquidus_dp = B2 + 2 * B3 * pressure;
      const double dT_liquidus_dp       = C2 + 2 * C3 * pressure;

      // We only consider peridotite melting (no pyroxenite)
      const double peridotite_fraction = 1.0;

      double melt_fraction_derivative_pressure = 0.;
      double melt_fraction_derivative_temperature = 0.;

      if (temperature > T_solidus && temperature < T_liquidus && pressure < 1.3e10)
        {
          // melt fraction when clinopyroxene is still present
          melt_fraction_derivative_temperature
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              / (T_lherz_liquidus - T_solidus);

          melt_fraction_derivative_pressure
            = beta * pow((temperature - T_solidus)/(T_lherz_liquidus - T_solidus),beta-1)
              * (dT_solidus_dp * (temperature - T_lherz_liquidus)
                 + dT_lherz_liquidus_dp * (T_solidus - temperature))
              / pow(T_lherz_liquidus - T_solidus,2);

          // melt fraction after melting of all clinopyroxene
          const double R_cpx = r1 + r2 * pressure;
          const double F_max = M_cpx / R_cpx;

          if (peridotite_melt_fraction(temperature, pressure, compositional_fields, position) > F_max)
            {
              const double T_max = std::pow(F_max,1.0/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
              const double dF_max_dp = - M_cpx * std::pow(r1 + r2 * pressure,-2) * r2;
              const double dT_max_dp = dT_solidus_dp
                                       + 1.0/beta * std::pow(F_max,1.0/beta - 1.0) * dF_max_dp * (T_lherz_liquidus - T_solidus)
                                       + std::pow(F_max,1.0/beta) * (dT_lherz_liquidus_dp - dT_solidus_dp);

              melt_fraction_derivative_temperature
                = (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  / (T_liquidus - T_max);

              melt_fraction_derivative_pressure
                = dF_max_dp
                  - dF_max_dp * std::pow((temperature - T_max)/(T_liquidus - T_max),beta)
                  + (1.0 - F_max) * beta * std::pow((temperature - T_max)/(T_liquidus - T_max),beta-1)
                  * (dT_max_dp * (T_max - T_liquidus) - (dT_liquidus_dp - dT_max_dp) * (temperature - T_max)) / std::pow(T_liquidus - T_max, 2);
            }

          double melt_fraction_derivative = 0;
          if (dependence == NonlinearDependence::temperature)
            melt_fraction_derivative = melt_fraction_derivative_temperature;
          else if (dependence == NonlinearDependence::pressure)
            melt_fraction_derivative = melt_fraction_derivative_pressure;
          else
            AssertThrow(false, ExcMessage("Error in calculating melt fraction derivative: not implemented"));

          entropy_gradient += melt_fraction_derivative * peridotite_melting_entropy_change * peridotite_fraction;
        }

      return entropy_gradient;
    }


    template <int dim>
    void
    PhaseTransitions<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &) const
    {
      for (unsigned int q=0; q<in.temperature.size(); ++q)
        melt_fraction(in.temperature[q],
                      std::max(0.0, in.pressure[q]),
                      in.composition[q],
                      in.position[q]);
      return;
    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    peridotite_melt_fraction (const double temperature,
                              const double pressure,
                              const std::vector<double> &,
                              const Point<dim> &) const
    {
      // anhydrous melting of peridotite after Katz, 2003
      const double T_solidus  = A1 + 273.15
                                + A2 * pressure
                                + A3 * pressure * pressure;
      const double T_lherz_liquidus = B1 + 273.15
                                      + B2 * pressure
                                      + B3 * pressure * pressure;
      const double T_liquidus = C1 + 273.15
                                + C2 * pressure
                                + C3 * pressure * pressure;

      // melt fraction for peridotite with clinopyroxene
      double peridotite_melt_fraction;
      if (temperature < T_solidus || pressure > 1.3e10)
        peridotite_melt_fraction = 0.0;
      else if (temperature > T_lherz_liquidus)
        peridotite_melt_fraction = 1.0;
      else
        peridotite_melt_fraction = std::pow((temperature - T_solidus) / (T_lherz_liquidus - T_solidus),beta);

      // melt fraction after melting of all clinopyroxene
      const double R_cpx = r1 + r2 * pressure;
      const double F_max = M_cpx / R_cpx;

      if (peridotite_melt_fraction > F_max && temperature < T_liquidus)
        {
          const double T_max = std::pow(F_max,1/beta) * (T_lherz_liquidus - T_solidus) + T_solidus;
          peridotite_melt_fraction = F_max + (1 - F_max) * pow((temperature - T_max) / (T_liquidus - T_max),beta);
        }
      return peridotite_melt_fraction;

    }

    template <int dim>
    double
    PhaseTransitions<dim>::
    melt_fraction (const double temperature,
                   const double pressure,
                   const std::vector<double> &composition, /*composition*/
                   const Point<dim> &position) const
    {
      return peridotite_melt_fraction(temperature, pressure, composition, position);

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
            prm.declare_entry ("Constant grain size", "0.003",Patterns::List (Patterns::Double(0)),
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
            prm.declare_entry ("Use dislocation creep", "false",
                               Patterns::Bool (),
                               "Whether to calculate viscosity with dislocation creep or not ");
            prm.declare_entry ("Use temperature jumps", "true",
                               Patterns::Bool (),
                               "Whether to calculate viscosity with dislocation creep or not ");


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
            prm.declare_entry ("Phase tracker", "0",
                               Patterns::List (Patterns::Double(0)),
                               "Used to track phase with compositional fields.");

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

            // Melting parameters
            prm.declare_entry ("Thermal expansion coefficient of melt", "6.8e-5",
                               Patterns::Double (0),
                               "The value of the thermal expansion coefficient $\\alpha_f$. "
                               "Units: $1/K$.");
          prm.declare_entry ("A1", "1085.7",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the solidus "
                             "of peridotite. "
                             "Units: $°C$.");
          prm.declare_entry ("A2", "1.329e-7",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("A3", "-5.1e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the solidus of peridotite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("B1", "1475.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the lherzolite "
                             "liquidus used for calculating the fraction "
                             "of peridotite-derived melt. "
                             "Units: $°C$.");
          prm.declare_entry ("B2", "8.0e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("B3", "-3.2e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the  lherzolite liquidus used for "
                             "calculating the fraction of peridotite-"
                             "derived melt. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("C1", "1780.0",
                             Patterns::Double (),
                             "Constant parameter in the quadratic "
                             "function that approximates the liquidus "
                             "of peridotite. "
                             "Units: $°C$.");
          prm.declare_entry ("C2", "4.50e-8",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $°C/Pa$.");
          prm.declare_entry ("C3", "-2.0e-18",
                             Patterns::Double (),
                             "Prefactor of the quadratic pressure term "
                             "in the quadratic function that approximates "
                             "the liquidus of peridotite. "
                             "Units: $°C/(Pa^2)$.");
          prm.declare_entry ("r1", "0.5",
                             Patterns::Double (),
                             "Constant in the linear function that "
                             "approximates the clinopyroxene reaction "
                             "coefficient. "
                             "Units: non-dimensional.");
          prm.declare_entry ("r2", "8e-11",
                             Patterns::Double (),
                             "Prefactor of the linear pressure term "
                             "in the linear function that approximates "
                             "the clinopyroxene reaction coefficient. "
                             "Units: $1/Pa$.");
          prm.declare_entry ("beta", "1.5",
                             Patterns::Double (),
                             "Exponent of the melting temperature in "
                             "the melt fraction calculation. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Peridotite melting entropy change", "-300",
                             Patterns::Double (),
                             "The entropy change for the phase transition "
                             "from solid to melt of peridotite. "
                             "Units: $J/(kg K)$.");
          prm.declare_entry ("Mass fraction cpx", "0.15",
                             Patterns::Double (),
                             "Mass fraction of clinopyroxene in the "
                             "peridotite to be molten. "
                             "Units: non-dimensional.");
          prm.declare_entry ("Relative density of melt", "0.9",
                             Patterns::Double (),
                             "The relative density of melt compared to the "
                             "solid material. This means, the density change "
                             "upon melting is this parameter times the density "
                             "of solid material."
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
          use_T_jumps                   = prm.get_bool ("Use temperature jumps");

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
          advect_log_grainsize                   = prm.get_bool ("Advect logarithm of grain size");
          recrystallized_grain_size = Utilities::string_to_double
                                      (Utilities::split_string_list(prm.get ("Recrystallized grain size")));

          //viscosity parameters
          dislocation_viscosity_iteration_threshold = prm.get_double("Dislocation viscosity iteration threshold");
          dislocation_viscosity_iteration_number = prm.get_integer("Dislocation viscosity iteration number");
          max_temperature_dependence_of_eta     = prm.get_double ("Maximum temperature dependence of viscosity");
          max_eta                              = prm.get_double ("Maximum viscosity");
          min_eta                              = prm.get_double("Minimum viscosity");
          constant_grain_size                  = Utilities::string_to_double
                                                 (Utilities::split_string_list(prm.get ("Constant grain size")));
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
          use_dislocation                   = prm.get_bool ("Use dislocation creep");

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
          phase_track                = Utilities::string_to_double
                                         (Utilities::split_string_list(prm.get ("Phase tracker")));

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

          // Melting section
          melt_thermal_alpha         = prm.get_double ("Thermal expansion coefficient of melt");
          A1              = prm.get_double ("A1");
          A2              = prm.get_double ("A2");
          A3              = prm.get_double ("A3");
          B1              = prm.get_double ("B1");
          B2              = prm.get_double ("B2");
          B3              = prm.get_double ("B3");
          C1              = prm.get_double ("C1");
          C2              = prm.get_double ("C2");
          C3              = prm.get_double ("C3");
          r1              = prm.get_double ("r1");
          r2              = prm.get_double ("r2");
          beta            = prm.get_double ("beta");
          peridotite_melting_entropy_change
            = prm.get_double ("Peridotite melting entropy change");
          relative_melt_density = prm.get_double ("Relative density of melt");
          M_cpx           = prm.get_double ("Mass fraction cpx");

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

    template <int dim>
    void
    PhaseTransitions<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<PhaseAdditionalOutputs<dim> >() == NULL)
        {
          const unsigned int n_points = out.viscosities.size();
          out.additional_outputs.push_back(
            std::shared_ptr<MaterialModel::AdditionalMaterialOutputs<dim> >
            (new MaterialModel::PhaseAdditionalOutputs<dim> (n_points)));
        }
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

