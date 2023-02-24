/*
  Copyright (C) 2019 - 2022 by the authors of the ASPECT code.

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

#include <aspect/material_model/rheology/pore_fluid_pressure.h>
#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <aspect/utilities.h>
#include <aspect/postprocess/particles.h>
#include <aspect/particle/property/interface.h>
#include <aspect/simulator.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <aspect/material_model/utilities.h>


namespace aspect
{
  namespace MaterialModel
  {
    namespace Rheology
    {

      template <int dim>
      double
      PoreFluidPressure<dim>::
      compute_fluid_ratio(const unsigned int composition,
                          const double depth) const
      {
        // Create a linear function to determine fluid ratio decrease with depth
        double current_fluid_ratio = 0;

        //if(this->get_timestep_number() > 0)
        //  std::cout<<depth<<std::endl;

        if(depth <= fluid_ratio_length[composition])
          current_fluid_ratio = fluid_ratio_top[composition] + depth*(fluid_ratio_base[composition] - fluid_ratio_top[composition])/fluid_ratio_length[composition];
        else if(depth <= fluid_ratio_cutoff)
          current_fluid_ratio = fluid_ratio_base[composition];
        else if(depth <= fluid_ratio_cutoff + fluid_cutoff_taper)
          current_fluid_ratio = fluid_ratio_base[composition] - (depth - fluid_ratio_cutoff)*(fluid_ratio_base[composition])/(fluid_cutoff_taper);


        return current_fluid_ratio;
      }


      template <int dim>
      void
      PoreFluidPressure<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Include pore fluid pressure", "false",
                           Patterns::Bool (),
                           "Whether to include Peierls creep in the rheological formulation.");
        prm.declare_entry ("Fluid ratio top", "0", Patterns::Anything(),
                           "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}. "
                           "List with as many components as active "
                           "compositional fields (material data is assumed to "
                           "be in order with the ordering of the fields). ");
        prm.declare_entry ("Fluid ratio base", "0", Patterns::Anything(),
                           "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}. "
                           "List with as many components as active "
                           "compositional fields (material data is assumed to "
                           "be in order with the ordering of the fields). ");
        prm.declare_entry ("Fluid ratio linear depth", "0", Patterns::Anything(),
                           "Lower cutoff for effective viscosity. Units: \\si{\\pascal\\second}. "
                           "List with as many components as active "
                           "compositional fields (material data is assumed to "
                           "be in order with the ordering of the fields). ");
        prm.declare_entry ("Pore fluid cutoff depth", "100e3", Patterns::Double (0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
        prm.declare_entry ("Pore fluid cutoff taper", "25e3", Patterns::Double (0.),
                           "Stabilizes strain dependent viscosity. Units: \\si{\\per\\second}.");
      }

      template <int dim>
      void
      PoreFluidPressure<dim>::parse_parameters (ParameterHandler &prm,
                                                const std::unique_ptr<std::vector<unsigned int>> &expected_n_phases_per_composition)
      {

        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();
        // Establish that a background field is required here
        const bool has_background_field = true;

        fluid_ratio_top = Utilities::parse_map_to_double_array (prm.get("Fluid ratio top"),
                                                          list_of_composition_names,
                                                          has_background_field,
                                                          "Fluid ratio top",
                                                          true,
                                                          0);
        fluid_ratio_base = Utilities::parse_map_to_double_array (prm.get("Fluid ratio base"),
                                                          list_of_composition_names,
                                                          has_background_field,
                                                          "Fluid ratio base",
                                                          true,
                                                          0);
        fluid_ratio_length = Utilities::parse_map_to_double_array (prm.get("Fluid ratio linear depth"),
                                                          list_of_composition_names,
                                                          has_background_field,
                                                          "Fluid ratio linear depth",
                                                          true,
                                                          0);
        fluid_ratio_cutoff = prm.get_double("Pore fluid cutoff depth");
        fluid_cutoff_taper = prm.get_double("Pore fluid cutoff taper");
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
#define INSTANTIATE(dim) \
  namespace Rheology \
  { \
    template class PoreFluidPressure<dim>; \
  }

    ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
  }
}
