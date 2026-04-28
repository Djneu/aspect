/*
  Copyright (C) 2017 - 2023 by the authors of the ASPECT code.

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


#include "initial_melt_components.h"
#include <aspect/initial_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/melt.h>
#include <aspect/initial_composition/function.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    InitialMeltComponents<dim>::initialize()
    {
      // Make sure we keep track of the initial temperature manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_temperature_manager = this->get_initial_temperature_manager_pointer();

      // Make sure we keep track of the initial composition manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_composition_manager = this->get_initial_composition_manager_pointer();
    }



    template <int dim>
    double
    InitialMeltComponents<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {
      AssertThrow(this->introspection().compositional_name_exists("porosity"),
                  ExcMessage("The initial composition plugin `porosity' did not find a "
                             "compositional field called `porosity' to initialize. Please add a "
                             "compositional field with this name."));

      const unsigned int porosity_index = this->introspection().compositional_index_for_name("porosity");
      const unsigned int mcs_index = this->introspection().compositional_index_for_name("morb_cs");
      const unsigned int mcl_index = this->introspection().compositional_index_for_name("morb_cl");
      const unsigned int ccs_index = this->introspection().compositional_index_for_name("cmorb_cs");
      const unsigned int ccl_index = this->introspection().compositional_index_for_name("cmorb_cl");
      const unsigned int hcs_index = this->introspection().compositional_index_for_name("hmorb_cs");
      const unsigned int hcl_index = this->introspection().compositional_index_for_name("hmorb_cl");

        MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
        MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());
        in.requested_properties = MaterialModel::MaterialProperties::density;

        in.position[0] = position;
        in.temperature[0] = initial_temperature_manager->initial_temperature(position);
        in.pressure[0] = this->get_adiabatic_conditions().pressure(position);
        const double density = this->get_adiabatic_conditions().density(position);
        in.pressure_gradient[0] = 0.0;
        in.velocity[0] = 0.0;    
                  
        const Utilities::NaturalCoordinate<dim> point =
          this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);

        std::vector<double> composition(this->n_compositional_fields());
        for (unsigned int i = 0; i < this->n_compositional_fields(); ++i)
          composition[i] = function->value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),i);
        
        auto [vfrac, melt_reaction_rate, solids, liquids, enthalpy] 
            = volatile_model.equilibrium(composition, in.temperature[0], 
                                in.pressure[0], density, 0, 1);
        
        if (compositional_index == porosity_index) return vfrac;
        if (compositional_index == mcs_index) return solids[1];
        if (compositional_index == ccs_index) return solids[2];
        if (compositional_index == hcs_index) return solids[3];
        //if (compositional_index == mcl_index) return liquids[1];
        if (compositional_index == mcl_index)
            return (vfrac > 0.0 ? liquids[1] : 0.0);
        //if (compositional_index == ccl_index) return liquids[2];
        if (compositional_index == ccl_index)
            return (vfrac > 0.0 ? liquids[2] : 0.0);
        //if (compositional_index == hcl_index) return liquids[3];
        if (compositional_index == hcl_index)
            return (vfrac > 0.0 ? liquids[3] : 0.0);
    
      return std::numeric_limits<double>::quiet_NaN();
    }

    template <int dim>
    void
    InitialMeltComponents<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        volatile_model.initialize_simulator (this->get_simulator());
        volatile_model.parse_parameters(prm);
        prm.leave_subsection();
      }
      prm.leave_subsection();


      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Function");
        {
          coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        }

        try
          {
            function
              = std::make_unique<Functions::ParsedFunction<dim>>(this->n_compositional_fields());
            function->parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Initial composition model.Function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'\n"
                      << "More information about the cause of the parse error \n"
                      << "is shown below.\n";
            throw;
          }

        prm.leave_subsection();
      }
      prm.leave_subsection();   
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(InitialMeltComponents,
                                              "initial melt components",
                                              "A class that implements initial conditions for the porosity field "
                                              "by computing the equilibrium melt fraction for the given initial "
                                              "condition and reference pressure profile. Note that this plugin only "
                                              "works if there is a compositional field called `porosity', and the "
                                              "used material model implements the 'MeltFractionModel' interface. "
                                              "For all compositional fields except porosity this plugin returns 0.0, "
                                              "and they are therefore not changed as long as the default `add' "
                                              "operator is selected for this plugin.")
  }
}
