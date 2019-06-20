/*
  Copyright (C) 2018 by the authors of the ASPECT code.
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


#include <aspect/mesh_deformation/fastscape.h>
#include <aspect/geometry_model/box.h>
#include <deal.II/grid/grid_generator.h>

#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {
  template <int dim>
  void
  FastScape<dim>::initialize ()
  {
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      AssertThrow(geometry != nullptr,
                        ExcMessage("Fastscape can only be run with a box model"));
      AssertThrow(dim == 3,
                        ExcMessage("Fastscape can only be applied to 3D runs."));

      x_extent = geometry->get_extents()[0];
      y_extent = geometry->get_extents()[1];

      const std::pair<int, int > repetitions = geometry->get_repetitions();
      int x_repetitions = repetitions.first;
      int y_repetitions = repetitions.second;

      nx = 1+std::pow(2,refinement)*x_repetitions;
      ny = 1+std::pow(2,refinement)*y_repetitions;
      array_size = nx*ny-1;
  }


    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                                    ConstraintMatrix &mesh_velocity_constraints,
                                                                    const std::set<types::boundary_id> &boundary_ids) const
    {
          int istep;
          int steps = nstep;

          /*
           * Pointers get allocated only 8 bytes at compile so don't run into
           * allocation issues. This creates an array that can delete itself.
           */
          std::unique_ptr<double[]> h (new double[array_size]);
          std::unique_ptr<double[]> h_temp (new double[array_size]);
          std::unique_ptr<double[]> kf (new double[array_size]);
          std::unique_ptr<double[]> kd (new double[array_size]);

          /* Initialize kf and kd across array.
           * TODO: Find a cleaner way to initialize these to the same value,
           * and find how to initialize if they're spatially dependent.
           */
          for (int i=0; i<=array_size; i++)
            {
              kf[i]=kff;
              kd[i]=kdd;
              h_temp[i]=0;
            }

          const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");

          // Get a quadrature rule that exists only on the corners
          QTrapez<dim-1> face_corners;
          FEFaceValues<dim> face_vals (this->get_mapping(), this->get_fe(), face_corners, update_quadrature_points | update_JxW_values);

          /* loop over all of the surface cells and save the position and values
           * TODO: Get velocity components as well
           */
          typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = this->get_triangulation().begin_active(),
                                                                                   endc = this->get_triangulation().end();
          for (; cell != endc; ++cell)
            if (cell->is_locally_owned() && cell->at_boundary())
              for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                if (cell->face(face_no)->at_boundary())
                  {
                    if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                      continue;

                    face_vals.reinit( cell, face_no);

                    //Find out dx and dy, and the total number of x points per row.
                    double dx = cell->extent_in_direction(0);
                    double dy = cell->extent_in_direction(1);
                    double numx = 1+x_extent/dx;

                    for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                      {
                        const Point<dim> vertex = face_vals.quadrature_point(corner);
                        const double elevation = vertex(dim-1);

                        double indx = 1+vertex(0)/dx;
                        double indy = 1+vertex(1)/dy;
                        double index = (indy-1)*numx+indx;

                        h_temp[index-1] = elevation;
                      }
                  }

          double a_dt = this->get_timestep();

          double f_dt = a_dt/steps;
          while (f_dt>max_timestep)
          {
        	  steps=steps*2;
        	  f_dt = a_dt/steps;
          }


          if(Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
          {
        	  //Set h to h values from processor zero
              for (int i=0; i<=array_size; i++)
            	  h[i]=h_temp[i];

              for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
                {
                  MPI_Status status;
                  MPI_Recv(h_temp.get(), array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

                  /*
                   * To keep indexing correct, each proc fills h_temp with zeroes in locations that
                   * it doesn't own. This checks whether or not we have an actual value.
                   */
                  for (int i=0; i<=array_size; i++)
                	  if(h_temp[i]>0)
                		  h[i]=h_temp[i];

                }

            //Initialize fastscape
            fastscape_init_();
            fastscape_set_nx_ny_(&nx,&ny);
            fastscape_setup_();

            //set x and y extent
            fastscape_set_xl_yl_(&x_extent,&y_extent);

            //Set time step
            fastscape_set_dt_(&f_dt);

            //Initialize topography
            fastscape_init_h_(h.get());

            //Set erosional parameters
            fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);

            //set boundary conditions
            fastscape_set_bc_(&bc);

            //Initialize first time step
            fastscape_get_step_(&istep);
            fastscape_vtk_(h.get(), &vexp);

           /* if(a_dt > 0)
            {
            do
              {
                //execute step, this increases timestep counter
                fastscape_execute_step_();
                //get value of time step counter
                fastscape_get_step_(&istep);
                //outputs new h values
                //fastscape_copy_h_(h.get());
                //output vtk
                fastscape_vtk_(h.get(), &vexp);
              }
            while (istep<steps);
            }*/

            //output timing
            fastscape_debug_();

            //end FastScape run
            fastscape_destroy_();
          }
          else
            {
              MPI_Send(h_temp.get(), array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());
            }
    }



    template <int dim>
    void FastScape<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscape");
        {
          prm.declare_entry("Number of steps", "5",
                            Patterns::Integer(),
	                        "Number of steps per ASPECT timestep");
		  prm.declare_entry("Maximum timestep", "10e3",
					         Patterns::Double(0),
                            "Maximum timestep for fastscape.");
          prm.declare_entry("Boundary conditions", "1111",
                            Patterns::Integer(),
                            "Boundary conditions where 0 is reflective and 1 is fixed  "
                            "height. Must be given in four digits, where the order is bottom,"
                            "right, top, left.");
          prm.declare_entry("Vertical exaggeration", "3",
                            Patterns::Double(),
                            "Vertical exaggeration for fastscape's VTK file.");
          prm.declare_entry("Maximum aspect refinement", "1",
                            Patterns::Integer(),
                            "Refinement level expected at surface to determine"
                            "proper nx and ny values");

          prm.enter_subsection ("Erosional parameters");
          {
            prm.declare_entry("Drainage area exponent", "0.5",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Slope exponent", "1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Multi-direction slope exponent", "1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Bedrock deposition coefficient", "0",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Sediment deposition coefficient", "0",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Bedrock river incision rate", "2e-6",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Sediment river incision rate", "-1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Bedrock diffusivity", "1e-1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
            prm.declare_entry("Sediment diffusivity", "-1",
                              Patterns::Double(),
                              "Theta parameter described in \\cite{KMM2010}. "
                              "An unstabilized free surface can overshoot its ");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();
    }

    template <int dim>
    void FastScape<dim>::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscape");
        {
          nstep = prm.get_integer("Number of steps");
          max_timestep = prm.get_double("Maximum timestep");
          bc = prm.get_integer("Boundary conditions");
          vexp = prm.get_double("Vertical exaggeration");
          nstep = prm.get_integer("Maximum aspect refinement");
          refinement = prm.get_integer("Maximum aspect refinement");

          prm.enter_subsection("Erosional parameters");
          {
            m = prm.get_double("Drainage area exponent");
            n = prm.get_double("Slope exponent");
            kfsed = prm.get_double("Sediment river incision rate");
            kff = prm.get_double("Bedrock river incision rate");
            kdsed = prm.get_double("Sediment diffusivity");
            kdd = prm.get_double("Bedrock diffusivity");
            g = prm.get_double("Bedrock deposition coefficient");
            gsed = prm.get_double("Sediment deposition coefficient");
            p = prm.get_double("Multi-direction slope exponent");
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();

      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation of the functions we implement in this file
namespace aspect
{
  namespace MeshDeformation
  {
    ASPECT_REGISTER_MESH_DEFORMATION_MODEL(FastScape,
                                           "fastscape",
                                           "A plugin, which prescribes the surface mesh to "
                                           "deform according to an analytically prescribed "
                                           "function. Note that the function prescribes a "
                                           "deformation velocity, i.e. the return value of "
                                           "this plugin is later multiplied by the time step length "
                                           "to compute the displacement increment in this time step. "
                                           "The format of the "
                                           "functions follows the syntax understood by the "
                                           "muparser library, see Section~\\ref{sec:muparser-format}.")
  }
}
