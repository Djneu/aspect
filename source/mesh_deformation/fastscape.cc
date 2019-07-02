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
#include <aspect/utilities.h>
#include <deal.II/base/table.h>
#include <deal.II/base/table_indices.h>

#include <deal.II/numerics/vector_tools.h>


#include <aspect/simulator.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/postprocess/particles.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>

namespace aspect
{
  namespace MeshDeformation
  {
    /*template <int dim>
     FastScape<dim>::FastScape()
       :
       function(dim)
     {}*/


    template <int dim>
    void
    FastScape<dim>::initialize ()
    {
      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      AssertThrow(geometry != nullptr,
                  ExcMessage("Fastscape can only be run with a box model"));
      AssertThrow((dy_slices & 1) != 0,
                 ExcMessage("Number of slices in y must be an odd number."));

      //second is for maximum coordiantes, first for minimum.f
      for (unsigned int i=0; i<dim; ++i)
      {
       grid_extent[i].first = 0;
       grid_extent[i].second = geometry->get_extents()[i];
      }


      //TODO: There has to be a better type to use to get this.
      const std::pair<int, int > repetitions = geometry->get_repetitions();
      int x_repetitions = repetitions.first;
      int y_repetitions = repetitions.second;

      /*
       * Set nx and ny based on repetitions and initial global refinement level in model.
       * NOTE: fastscape refinement must be at least as high as the refinement on surface.
       * TODO: There needs to be a good error message if refinement is chosen
       * incorrectly. Right now it just throws a segmentation fault. Also,
       * should I add adapative here? This needs to actually equal the refinement
       * level at the surface.
       */
      nx = 1+std::pow(2,initial_global_refinement+additional_refinement)*x_repetitions;
      dx = grid_extent[0].second/(nx-1);

      //sub intervals are 1 less than points.
      table_intervals[0] = nx-1;
      //TODO: it'd be best to not have to use dim-1 intervals at all.
      table_intervals[dim-1] = 1;

  if (dim == 2)
  {
      ny = dy_slices;
      dy = dx;
      grid_extent[1].second = (ny-1)*dy;

  }

  if (dim == 3)
  {
      ny = 1+std::pow(2,initial_global_refinement+additional_refinement)*y_repetitions;
      dy = grid_extent[1].second/(ny-1);
      table_intervals[1] = ny-1;
  }


      array_size = nx*ny-1;
      numx = 1+grid_extent[0].second/dx;
    }


    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                             ConstraintMatrix &mesh_velocity_constraints,
                                                             const std::set<types::boundary_id> &boundary_ids) const
    {
      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      const bool top_boundary = boundary_ids.find(relevant_boundary) == boundary_ids.begin();
      double a_dt = this->get_old_timestep()/year_in_seconds;


      if (a_dt > 0 && top_boundary)
        {

          /*
           * Initialize a vector of temporary variables to hold: z component, Vx, Vy, and Vz.
           * For some reason I need to add one to the temporary array_size but not V.
           */
          std::vector<std::vector<double>> temporary_variables(dim+1, std::vector<double>(array_size+1));
          std::vector<double> V(array_size);


          // Get a quadrature rule that exists only on the corners

          const QIterated<dim-1> face_corners (QTrapez<1>(),
                                                   pow(2,additional_refinement));

          FEFaceValues<dim> fe_face_values (this->get_mapping(),
                                            this->get_fe(),
                                            face_corners,
                                            update_values |
                                            update_quadrature_points);

          typename DoFHandler<dim>::active_cell_iterator
          cell = this->get_dof_handler().begin_active(),
          endc = this->get_dof_handler().end();

          for (; cell != endc; ++cell)
            if (cell->is_locally_owned() && cell->at_boundary())
              for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
                if (cell->face(face_no)->at_boundary())
                  {
                    if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                      continue;

                    std::vector<Tensor<1,dim> > vel( face_corners.size() );
                    fe_face_values.reinit( cell, face_no);
                    fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel );

                    for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                      {
                        const Point<dim> vertex = fe_face_values.quadrature_point(corner);
                        double indx = 1+vertex(0)/dx;

                      if(dim == 2)
                      {
                          for (unsigned int ys=0; ys<dy_slices; ys++)
                          {
                           double index = indx+numx*ys;
                           temporary_variables[0][index-1] = this->get_geometry_model().height_above_reference_surface(vertex); //vertex(dim-1);   //z component

                           for (unsigned int i=0; i<dim; ++i)
                               temporary_variables[i+1][index-1] = vel[corner][i];
                          }
                      }

                      if(dim == 3)
                      {
                          double indy = 1+vertex(1)/dy;
                          double index = (indy-1)*numx+indx;
                          temporary_variables[0][index-1] = this->get_geometry_model().height_above_reference_surface(vertex); //vertex(dim-1);   //z component

                          for (unsigned int i=0; i<dim; ++i)
                              temporary_variables[i+1][index-1] = vel[corner][i];
                      }


                      }
                  }

          //Only run fastscape on a single processor
          if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
            {
              /* Initialize the variables that will be sent to fastscape.
               * These have to be doubles of array_size, which C++ doesn't like,
               * so they're initialized this way.
               */
              std::unique_ptr<double[]> h (new double[array_size]);
              std::unique_ptr<double[]> vx (new double[array_size]);
              std::unique_ptr<double[]> vy (new double[array_size]);
              std::unique_ptr<double[]> vz (new double[array_size]);
              std::unique_ptr<double[]> kf (new double[array_size]);
              std::unique_ptr<double[]> kd (new double[array_size]);
              int istep;
              int steps = nstep;

              /* Initialize kf and kd across array, and set h values to what proc zero has.
               * TODO: Find a cleaner way to set these to the same value,
               */
              for (int i=0; i<=array_size; i++)
                {
                  h[i]= temporary_variables[0][i];
                  vx[i]=temporary_variables[1][i];
                  vz[i]=temporary_variables[dim][i];

                  if(dim == 2 )
                      vy[i]=0;

                  if(dim == 3)
                	  vy[i]=temporary_variables[2][i];

                  kf[i] = kff;
                  kd[i] = kdd;
                }

              for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
                {
                  MPI_Status status;
                  for (int i=0; i<dim+1; i++)
                     MPI_Recv(&temporary_variables[i][0], array_size+1, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

                  /*
                   * Each processor initialized everything to zero, so if it had any
                   * points at the surface a value will be added.
                   */
                 for (int i=0; i<=array_size; i++)
                	if(temporary_variables[0][i]>0)
                    {
                        h[i]= temporary_variables[0][i];
                        vx[i]=temporary_variables[1][i];
                        vz[i]=temporary_variables[dim][i];

                        if(dim == 2 )
                            vy[i]=0;

                        if(dim == 3)
                      	    vy[i]=temporary_variables[2][i];
                    }

                }

              //Reset h in temporary to find velocities later on.
              for (int i=0; i<=array_size; i++)
                {
                  temporary_variables[0][i] = h[i];
                }

              double f_dt = a_dt/steps;
              while (f_dt>max_timestep)
                {
                  steps=steps*2;
                  f_dt = a_dt/steps;
                }
              std::cout<<"Fastscape timestep: "<<f_dt<<"  "<<a_dt<<"  "<<steps<<std::endl;

              //Initialize fastscape
              fastscape_init_();
              fastscape_set_nx_ny_(&nx,&ny);
              fastscape_setup_();

              //set x and y extent
             fastscape_set_xl_yl_(&grid_extent[0].second,&grid_extent[1].second);

              //Set time step
              //f_dt = 2e4;
              fastscape_set_dt_(&f_dt);

              //Initialize topography
              fastscape_init_h_(h.get());

              //Set erosional parameters TODO: why is it g twice here?
              fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);

              //set boundary conditions
              fastscape_set_bc_(&bc);

              //Initialize first time step
              fastscape_get_step_(&istep);

              fastscape_vtk_(h.get(), &vexp);

             do
                {
                  //execute step, this increases timestep counter
                  fastscape_execute_step_();
                  //get value of time step counter
                  fastscape_get_step_(&istep);
                  //outputs new h values
                  fastscape_copy_h_(h.get());
                  //output vtk
                  fastscape_vtk_(h.get(), &vexp);
                }
              while (istep<steps);

              //output timing
              fastscape_debug_();

              //end FastScape run
             fastscape_destroy_();

              for (int i=0; i<=array_size; i++)
              {
                V[i] = (h[i] - temporary_variables[0][i])/a_dt;
              }

              MPI_Bcast(&V[0], array_size+1, MPI_DOUBLE, 0, this->get_mpi_communicator());

            }
          else
            {
        	  for (int i=0; i<dim+1; i++)
                MPI_Send(&temporary_variables[i][0], array_size+1, MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

        	  MPI_Bcast(&V[0], array_size+1, MPI_DOUBLE, 0, this->get_mpi_communicator());
            }

          TableIndices<dim> size_idx;
          for (unsigned int d=0; d<dim; ++d)
          {
            size_idx[d] = table_intervals[d]+1;
          }

          Table<dim,double> data_table;
          data_table.TableBase<dim,double>::reinit(size_idx);
          TableIndices<dim> idx;

          //Average 2d slices, this may depend on chosen boundary conditions.
          if (dim == 2)
          {
              std::vector<double> V2(numx);

        	  for (int i=0; i<numx; i++)
        	  {
              for (unsigned int ys=0; ys<dy_slices; ys++)
              {
                int index = i+numx*ys;
                V2[i] += V[index];
              }

              V2[i] = V2[i]/dy_slices;
        	 }
                  for (unsigned int i=0; i<data_table.size()[1]; ++i)
                    {
                      idx[1] = i;
                      for (unsigned int j=0; j<data_table.size()[0]; ++j)
                        {
                          idx[0] = j;
                          if(i == 1 )
                            data_table(idx) = V2[j]/year_in_seconds;
                          else
                        	data_table(idx)= 0;
                        }

                    }

          }

          if(dim == 3)
          {
          //Indexes through y and then x
          for (unsigned int k=0; k<data_table.size()[2]; ++k)
            {
        	  idx[2] = k;
              for (unsigned int i=0; i<data_table.size()[1]; ++i)
                {
                  idx[1] = i;
                  for (unsigned int j=0; j<data_table.size()[0]; ++j)
                    {
                      idx[0] = j;
                      if(k==1 )
                        data_table(idx) = V[nx*i+j];
                      else
                    	data_table(idx)= 0;
                    }

                }
            }
          }


          Functions::InterpolatedUniformGridData<dim> *velocities;
              velocities = new Functions::InterpolatedUniformGridData<dim> (grid_extent,
            		                                                        table_intervals,
                                                                            data_table);

              auto lambda = [&](const Point<dim> &p) -> double
              {
                return velocities->value(p);
              };

          VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
            lambda,
            dim-1,
            dim);

            VectorTools::interpolate_boundary_values (mesh_deformation_dof_handler,
                                                      *boundary_ids.begin(),
                                                      vector_function_object,
                                                      mesh_velocity_constraints);

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
          prm.declare_entry("Additional fastscape refinement", "0",
                            Patterns::Integer(),
                            "Refinement level expected at surface to determine"
                            "proper nx and ny values");
          prm.declare_entry("Ny slices for 2d", "11",
                            Patterns::Integer(),
                            "Number of steps per ASPECT timestep");

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

      prm.enter_subsection ("Mesh refinement");
      {
        prm.declare_entry ("Initial global refinement", "2",
                           Patterns::Integer (0),
                           "The number of global refinement steps performed on "
                           "the initial coarse mesh, before the problem is first "
                           "solved there.");
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
          additional_refinement = prm.get_integer("Additional fastscape refinement");
          dy_slices = prm.get_integer("Ny slices for 2d");

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

      prm.enter_subsection ("Mesh refinement");
      {
        initial_global_refinement    = prm.get_integer ("Initial global refinement");
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