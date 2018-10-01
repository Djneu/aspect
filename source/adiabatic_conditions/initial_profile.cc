/*
  Copyright (C) 2011 - 2016 by the authors of the ASPECT code.

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


#include <aspect/global.h>
#include <aspect/material_model/interface.h>
#include <aspect/material_model/phase_transitions.h>
#include <aspect/adiabatic_conditions/initial_profile.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/initial_composition/interface.h>
#include <deal.II/base/signaling_nan.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace AdiabaticConditions
  {
    template <int dim>
    InitialProfile<dim>::InitialProfile()
      :
      initialized(false),
      n_points(1000),
      temperatures(n_points, numbers::signaling_nan<double>()),
      pressures(n_points, numbers::signaling_nan<double>()),
      densities(n_points, numbers::signaling_nan<double>())
    {}



    template <int dim>
    void
    InitialProfile<dim>::initialize()
    {
      if (initialized)
        return;

      Assert (pressures.size() == n_points, ExcInternalError());
      Assert (temperatures.size() == n_points, ExcInternalError());

      delta_z = this->get_geometry_model().maximal_depth() / (n_points-1);

      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());

      // Constant properties on the reference profile
      in.strain_rate.resize(0); // we do not need the viscosity
      in.velocity[0] = Tensor <1,dim> ();

      // Check whether gravity is pointing up / out or down / in. In the normal case it should
      // point down / in and therefore gravity should be positive, leading to increasing
      // adiabatic pressures and temperatures with depth. In some cases it will point up / out
      // (e.g. for backward advection), in which case the pressures and temperatures should
      // decrease with depth and therefore gravity has to be negative in the following equations.
      const Tensor <1,dim> g = this->get_gravity_model().gravity_vector(this->get_geometry_model().representative_point(0));
      const Point<dim> point_surf = this->get_geometry_model().representative_point(0);
      const Point<dim> point_bot = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
      const int gravity_direction =  (g * (point_bot - point_surf) >= 0) ?
                                     1 :
                                     -1;

      // Initialize variable for density
      double density = reference_rho;
      double density1 = reference_rho;

      // now integrate downward using the explicit Euler method for simplicity
      //
      // note: p'(z) = rho(p,T) * |g|
      //       T'(z) = alpha |g| T / C_p
      for (unsigned int i=0; i<n_points; ++i)
        {
 
          // Define current depth
          //const double z = double(i)/double(n_points-1)*this->get_geometry_model().maximal_depth(); 
          const double z = double(i)/double(n_points-1); 
          const double max_depth = this->get_geometry_model().maximal_depth();

          if (i==0)
            {
              pressures[i] = this->get_surface_pressure();
              temperatures[i] = this->get_adiabatic_surface_temperature();
              density = reference_rho;
            }
          else
            {

              // Using surface (reference) value of thermal expansivity
              const double alpha = 2.189e-5; //out.thermal_expansion_coefficients[0];

              // Handle the case that cp is zero (happens in simple Stokes test problems like sol_cx). By setting
              // 1/cp = 0.0 we will have a constant temperature profile with depth.
              const double one_over_cp = (out.specific_heat[0]>0.0) ? 1.0/out.specific_heat[0] : 0.0;

              // get the magnitude of gravity. we assume
              // that gravity always points along the depth direction. this
              // may not strictly be true always but is likely a good enough
              // approximation here.
              const double gravity = gravity_direction * this->get_gravity_model().gravity_vector(in.position[0]).norm();

              // Calculate Dissipation number and Grunheissen parameter 
              const double grun = alpha / ( reference_compressibility * reference_rho * reference_specific_heat );
              const double di = ( alpha * gravity * this->get_geometry_model().maximal_depth() ) / reference_specific_heat; 

              // Calculate current pressure (using previous density)
              pressures[i] = pressures[i-1] + density1 * gravity * delta_z;

              double phase_function_rho = 0.;
              double phase_function_t = 0.;
              double entropy_change = 0.;
              double temperature_jump = 0.;

      
              for (unsigned int j=0; j<transition_depths.size(); ++j)
                {


                  double depth_deviation = (z*max_depth) - transition_depths[j];


                  phase_function_rho += 0.5 * ( 1. + tanh( depth_deviation / transition_widths[j] ) )*density_jumps[j];


                  entropy_change = transition_slopes[j] * density_jumps[j] / std::pow(reference_rho, 2);


                  temperature_jump = this->get_adiabatic_surface_temperature() * std::exp( di * transition_depths[j] / max_depth )
                                       * entropy_change * one_over_cp;
                  phase_function_t += 0.5 * ( 1. + tanh( depth_deviation / transition_widths[j] ) )*temperature_jump;      
                }


              // Recalculate density. Divide by "1" is Gruenheissen parameter
             // density = reference_rho * std::pow( ( ( (delta * di * z ) / 1. ) + 1. ), 1 / delta ) +  phase_function_rho;
              density1 = reference_rho * std::exp( di * z / grun ) + phase_function_rho;
              density = reference_rho * std::exp( di * z / grun );

              //calculate the temperature profile with or without temperature jumps associated with phase transitions.
              if(use_T_jumps)  
              {
              temperatures[i] = (this->include_adiabatic_heating())
                                ?
                                this->get_adiabatic_surface_temperature() * std::exp( di * z) +  phase_function_t
                                :
                                this->get_adiabatic_surface_temperature();
              }
              else
              temperatures[i] = (this->include_adiabatic_heating())
                                ?
                                this->get_adiabatic_surface_temperature() * std::exp( di * z)
                                :
                                this->get_adiabatic_surface_temperature();
            }

          //const double z = double(i)/double(n_points-1)*this->get_geometry_model().maximal_depth();
          const Point<dim> representative_point = this->get_geometry_model().representative_point (z);
          const Tensor <1,dim> g = this->get_gravity_model().gravity_vector(representative_point);

          in.position[0] = representative_point;
          in.temperature[0] = temperatures[i];
          in.pressure[0] = pressures[i];

          // we approximate the pressure gradient by extrapolating the values
          // from the two points above
          if (i>0)
            in.pressure_gradient[0] = g/(g.norm() != 0.0 ? g.norm() : 1.0)
                                      * (pressures[i] - pressures[i-1]) / delta_z;
          else
            in.pressure_gradient[0] = Tensor <1,dim> ();

          if (reference_composition == initial_composition)
            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              in.composition[0][c] = this->get_initial_composition_manager().initial_composition(representative_point, c);
          else if (reference_composition == reference_function)
            {
              const double depth = this->get_geometry_model().depth(representative_point);
              const Point<1> p(depth);
              for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
                in.composition[0][c] = composition_function->value(p, c);
            }
          else
            AssertThrow(false,ExcNotImplemented());

          this->get_material_model().evaluate(in, out);

          densities[i] = density;
        }

      if (gravity_direction == 1 && this->get_surface_pressure() >= 0)
        {
          Assert (*std::min_element (pressures.begin(), pressures.end()) >=
                  -std::numeric_limits<double>::epsilon() * pressures.size(),
                  ExcMessage("Adiabatic InitialProfile encountered a negative pressure of "
                             + dealii::Utilities::to_string(*std::min_element (pressures.begin(), pressures.end()))));
        }
      else if (gravity_direction == -1 && this->get_surface_pressure() <= 0)
        {
          Assert (*std::max_element (pressures.begin(), pressures.end()) <=
                  std::numeric_limits<double>::epsilon() * pressures.size(),
                  ExcMessage("Adiabatic InitialProfile encountered a positive pressure of "
                             + dealii::Utilities::to_string(*std::max_element (pressures.begin(), pressures.end()))));
        }

      Assert (*std::min_element (temperatures.begin(), temperatures.end()) >=
              -std::numeric_limits<double>::epsilon() * temperatures.size(),
              ExcMessage("Adiabatic InitialProfile encountered a negative temperature."));


      initialized = true;
    }



    template <int dim>
    bool
    InitialProfile<dim>::is_initialized() const
    {
      return initialized;
    }



    template <int dim>
    double InitialProfile<dim>::pressure (const Point<dim> &p) const
    {
      return get_property(p,pressures);
    }



    template <int dim>
    double InitialProfile<dim>::temperature (const Point<dim> &p) const
    {
      return get_property(p,temperatures);
    }



    template <int dim>
    double InitialProfile<dim>::density (const Point<dim> &p) const
    {
      return get_property(p,densities);
    }



    template <int dim>
    double InitialProfile<dim>::density_derivative (const Point<dim> &p) const
    {
      const double z = this->get_geometry_model().depth(p);

      if (z >= this->get_geometry_model().maximal_depth())
        {
          Assert (z <= this->get_geometry_model().maximal_depth() + delta_z,
                  ExcInternalError());
          return (densities.back() - densities[densities.size()-2]) / delta_z;
        }

      if (z < 0)
        {
          Assert (z >= -delta_z, ExcInternalError());
          return (densities[1] - densities.front()) / delta_z;
        }

      // if z/delta_z is within [k-eps, k+eps] of a whole number k, round it down to k-1
      const unsigned int i = static_cast<unsigned int>((z/delta_z) * (1. - 2. * std::numeric_limits<double>::epsilon()));
      Assert (i < densities.size() - 1, ExcInternalError());

      return (densities[i+1]-densities[i])/delta_z;
    }



    template <int dim>
    double InitialProfile<dim>::get_property (const Point<dim> &p,
                                              const std::vector<double> &property) const
    {
      const double z = this->get_geometry_model().depth(p);

      if (z >= this->get_geometry_model().maximal_depth())
        {
          Assert (z <= this->get_geometry_model().maximal_depth() + delta_z,
                  ExcInternalError());
          return property.back();
        }

      if (z <= 0)
        {
          Assert (z >= -delta_z, ExcInternalError());
          return property.front();
        }

      const double floating_index = z/delta_z;
      const unsigned int i = static_cast<unsigned int>(floating_index);

      // If p is close to an existing value use that one. This prevents
      // asking for values at i+1 while initializing i+1 (when p is at the
      // depth of index i).
      if (std::abs(floating_index-std::floor(floating_index+0.5)) < 1e-6)
        return property[i];

      Assert (i+1 < property.size(), ExcInternalError());

      // now do the linear interpolation
      const double d = floating_index - i;
      Assert ((d>=0) && (d<=1), ExcInternalError());

      return d*property[i+1] + (1.-d)*property[i];
    }




    template <int dim>
    void
    InitialProfile<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Initial profile");
        {
          prm.declare_entry("Delta", "2",
                            Patterns::Double (0),
                            "Power law exponent of thermal expansivity. Units: none$.");

          Functions::ParsedFunction<1>::declare_parameters (prm, 1);
          prm.declare_entry("Composition reference profile","initial composition",
                            Patterns::Selection("initial composition|function"),
                            "Select how the reference profile for composition "
                            "is computed. This profile is used to evaluate the "
                            "material model, when computing the pressure and "
                            "temperature profile.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    InitialProfile<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Adiabatic conditions model");
      {
        prm.enter_subsection("Initial profile");
        {

          delta = prm.get_double("Delta");
          const std::string composition_profile = prm.get("Composition reference profile");

          if (composition_profile == "initial composition")
            reference_composition = initial_composition;
          else if (composition_profile == "function")
            reference_composition = reference_function;
          else
            AssertThrow(false, ExcNotImplemented());

          if ((this->n_compositional_fields() > 0) && (reference_composition == reference_function))
            {
              composition_function.reset(new Functions::ParsedFunction<1>(this->n_compositional_fields()));
              try
                {
                  composition_function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Adiabatic conditions model.Initial profile'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      // Get data from Material model parameters
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Phase transitions");
        {
          reference_rho              = prm.get_double ("Reference density");
          reference_T                = prm.get_double ("Reference temperature");
          reference_compressibility  = prm.get_double ("Compressibility");
          reference_specific_heat    = prm.get_double ("Reference specific heat");
          thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
          use_T_jumps                   = prm.get_bool ("Use temperature jumps");
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
        }
        prm.leave_subsection();
        //ASSERT(prm.get(
      }
      prm.leave_subsection();
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace AdiabaticConditions
  {
    ASPECT_REGISTER_ADIABATIC_CONDITIONS_MODEL(InitialProfile,
                                               "initial profile",
                                               "A model in which the adiabatic profile is "
                                               "calculated once at the start of the model run. "
                                               "The gravity is assumed to be in depth direction "
                                               "and the composition is either given by the initial "
                                               "composition at reference points or computed "
                                               "as a reference depth-function. "
                                               "All material parameters are computed by the "
                                               "material model plugin.")
  }
}
