/*
  Copyright (C) 2022 by the authors of the ASPECT code.
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

#ifdef ASPECT_WITH_FASTSCAPE

#include <aspect/mesh_deformation/fastscape.h>
#include <aspect/geometry_model/box.h>
#include <deal.II/numerics/vector_tools.h>
#include <aspect/postprocess/visualization.h>
#include <ctime>

namespace aspect
{
  namespace MeshDeformation
  {

    template <int dim>
    void
    FastScape<dim>::initialize ()
    {
      AssertThrow(Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("FastScape can only be run with a box geometry model."));

      const GeometryModel::Box<dim> *geometry
        = dynamic_cast<const GeometryModel::Box<dim>*> (&this->get_geometry_model());

      // Find the id associated with the top boundary and boundaries that call mesh deformation.
      const types::boundary_id top_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      const std::set<types::boundary_id> mesh_deformation_boundary_ids
        = this->get_mesh_deformation_handler().get_active_mesh_deformation_boundary_indicators();

      // Get the deformation type names called for each boundary.
      std::map<types::boundary_id, std::vector<std::string>> mesh_deformation_boundary_indicators_map
        = this->get_mesh_deformation_handler().get_active_mesh_deformation_names();

      // Loop over each mesh deformation boundary, and make sure FastScape is only called on the surface.
      for (std::set<types::boundary_id>::const_iterator p = mesh_deformation_boundary_ids.begin();
           p != mesh_deformation_boundary_ids.end(); ++p)
        {
          const std::vector<std::string> names = mesh_deformation_boundary_indicators_map[*p];
          for (unsigned int i = 0; i < names.size(); ++i )
            {
              if (names[i] == "fastscape")
                AssertThrow((*p == top_boundary),
                            ExcMessage("FastScape can only be called on the surface boundary."));
            }
        }

      // Initialize parameters for restarting FastScape
      restart = this->get_parameters().resume_computation;
      restart_step = 0;

      // Since we don't open these until we're on one process, we need to check if the
      // restart files exist before hand.
      // TODO: This was quickly done and can likely be shortened/improved.
      if (restart)
        {
          // Create variables for output directory and restart file
          std::string dirname = this->get_output_directory();

          std::ifstream in;
          in.open(dirname + "fastscape_h_restart.txt");
          if (in.fail())
            AssertThrow(false,ExcMessage("Cannot open topography file to restart FastScape."));
          in.close();

          in.open(dirname + "fastscape_b_restart.txt");
          if (in.fail())
            AssertThrow(false,ExcMessage("Cannot open basement file to restart FastScape."));
          in.close();

          in.open(dirname + "fastscape_steps_restart.txt");
          if (in.fail())
            AssertThrow(false,ExcMessage("Cannot open steps file to restart FastScape."));
          in.close();
        }

      // The first entry represents the minimum coordinates of the model domain, the second the maximum coordinates.
      for (unsigned int i=0; i<dim; ++i)
        {
          grid_extent[i].first = geometry->get_origin()[i];
          grid_extent[i].second = geometry->get_extents()[i];
        }

      // Get the x and y repetitions used in the parameter file so
      // the FastScape cell size can be properly set.
      std::array<unsigned int,dim> repetitions = geometry->get_repetitions();

      // Set number of x points, which is generally 1+(FastScape refinement level)^2.
      // The FastScape refinement level is a combination of the maximum ASPECT refinement level
      // at the surface and any additional refinement we want in FastScape. If
      // repetitions are specified we need to adjust the number of points to match what ASPECT has,
      // which can be determined by multiplying the points by the repetitions before adding 1.
      // Finally, if ghost nodes are used we add two additional points on each side.
      nx = 1+(2*use_ghost_nodes)+std::pow(2,maximum_surface_refinement_level+additional_refinement)*repetitions[0];

      // Size of FastScape cell.
      dx = (grid_extent[0].second)/(nx-1-(2*use_ghost_nodes));

      // FastScape X extent, which is generally ASPECT's extent unless the ghost nodes are used,
      // in which case 2 cells are added on either side.
      x_extent = (grid_extent[0].second)+2*dx*use_ghost_nodes;

      // Sub intervals are 3 less than points, if including the ghost nodes. Otherwise 1 less.
      table_intervals[0] = nx-1-(2*use_ghost_nodes);
      table_intervals[dim-1] = 1;

      if (dim == 2)
        {
          dy = dx;
          y_extent = round(y_extent_2d/dy)*dy+2*dy*use_ghost_nodes;
          ny = 1+y_extent/dy;
        }
      else
        {
          ny = 1+(2*use_ghost_nodes)+std::pow(2,maximum_surface_refinement_level+additional_refinement)*repetitions[1];
          dy = (grid_extent[1].second)/(ny-1-(2*use_ghost_nodes));
          table_intervals[1] = ny-1-(2*use_ghost_nodes);
          y_extent = (grid_extent[1].second)+2*dy*use_ghost_nodes;
        }

      // Determine array size to send to FastScape
      array_size = nx*ny;

      // Create a folder for the FastScape visualization files.
      Utilities::create_directory (this->get_output_directory() + "VTK/",
                                   this->get_mpi_communicator(),
                                   false);

      last_output_time = 0;
    }


    template <int dim>
    void
    FastScape<dim>::compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                                             AffineConstraints<double> &mesh_velocity_constraints,
                                                             const std::set<types::boundary_id> &boundary_ids) const
    {
      if (this->get_timestep_number() == 0)
        return;

      TimerOutput::Scope timer_section(this->get_computing_timer(), "FastScape plugin");
      const int current_timestep = this->get_timestep_number ();
      const double a_dt = this->get_timestep()/year_in_seconds;

      //Vector to hold the velocities that represent the change to the surface.
      std::vector<double> V(array_size);

      // FastScape requires multiple specially defined and ordered variables sent to its functions. To make
      // the transfer of these down to one process easier, we first fill out a vector of temporary variables,
      // then when we get down to one process we use these temporary variables to fill the double arrays
      // in the order needed for FastScape.
      std::vector<std::vector<double>> temporary_variables = get_aspect_values();

      // Run FastScape on single process.
      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Initialize the variables that will be sent to FastScape.
          // These have to be doubles of array_size, which C++ doesn't like,
          // so they're initialized this way.
          std::unique_ptr<double[]> h (new double[array_size]());
          std::unique_ptr<double[]> vx (new double[array_size]());
          std::unique_ptr<double[]> vy (new double[array_size]());
          std::unique_ptr<double[]> vz (new double[array_size]());
          std::unique_ptr<double[]> kf (new double[array_size]());
          std::unique_ptr<double[]> kd (new double[array_size]());
          std::unique_ptr<double[]> b (new double[array_size]());
          std::vector<double> h_old(array_size);

          // Create variables for output directory and restart file
          std::string dirname = this->get_output_directory();
          const std::string restart_filename = dirname + "fastscape_h_restart.txt";
          const std::string restart_step_filename = dirname + "fastscape_steps_restart.txt";
          const std::string restart_filename_basement = dirname + "fastscape_b_restart.txt";

          fill_fastscape_arrays(h.get(), kd.get(), kf.get(), vx.get(), vy.get(), vz.get(), temporary_variables);

          if (current_timestep == 1 || restart)
            {
              this->get_pcout() << "   Initializing FastScape... " << (1+maximum_surface_refinement_level+additional_refinement) <<
                                " levels, cell size: " << dx << " m." << std::endl;

              // If we are restarting from a checkpoint, load h values for FastScape instead of using the ASPECT values.
              if (restart)
                {
                  this->get_pcout() << "      Loading FastScape restart file... " << std::endl;

                  // Load in h values.
                  std::ifstream in;
                  in.open(restart_filename.c_str());
                  if (in)
                    {
                      int line = 0;

                      while (line < array_size)
                        {
                          in >> h[line];
                          line++;
                        }

                      in.close();
                    }

                  // Load in b values.
                  std::ifstream in_b;
                  in_b.open(restart_filename_basement.c_str());
                  if (in_b)
                    {
                      int line = 0;

                      while (line < array_size)
                        {
                          in_b >> b[line];
                          line++;
                        }

                      in_b.close();
                    }

                  // Now load the FastScape istep at time of restart.
                  // Reinitializing FastScape always resets this to 0, so here
                  // we keep it in a separate variable to keep track for visualization files.
                  std::ifstream in_step;
                  in_step.open(restart_step_filename.c_str());
                  if (in_step)
                    {
                      in_step >> restart_step;
                      in_step.close();
                    }

                  restart = false;
                }

              initialize_fastscape(h.get(), b.get(), kd.get(), kf.get());
            }
          else
            {
              // If it isn't the first timestep we ignore initialization and instead copy all height values from FastScape.
              if (use_velocities)
                fastscape_copy_h_(h.get());
            }

          // Find the appropriate sediment rain based off the time interval.
          double time = this->get_time()/year_in_seconds;
          double sediment_rain = sediment_rain_rates[0];
          for (unsigned int j=0; j<sediment_rain_times.size(); j++)
            {
              if (time > sediment_rain_times[j])
                sediment_rain = sediment_rain_rates[j+1];
            }

          // Keep initial h values so we can calculate velocity later.
          // In the first timestep, h will be given from other processes.
          // In later timesteps, we copy h directly from FastScape.
          std::srand(fs_seed);
          for (int i=0; i<array_size; i++)
            {
              h_old[i] = h[i];

              // Initialize random noise after h_old is set, so aspect sees this initial topography change.
              if (current_timestep == 1)
                {
                  // + or - topography based on the initial noise magnitude.
                  const double h_seed = (std::rand()%( 2*noise_h+1 )) - noise_h;
                  h[i] = h[i] + h_seed;
                }

              // Here we add the sediment rain (m/yr) as a flat increase in height.
              // This is done because adding it as an uplift rate would affect the basement.
              if (sediment_rain > 0 && use_marine)
                {
                  // Only apply sediment rain to areas below sea level.
                  if (h[i] < sl)
                    {
                      // If the rain would put us above sea level, set height to sea level.
                      if (h[i] + sediment_rain*a_dt > sl)
                        h[i] = sl;
                      else
                        h[i] = std::min(sl,h[i] + sediment_rain*a_dt);
                    }
                }
            }

          // The ghost nodes are added as a single layer of points surrounding the entire model.
          // For example, if ASPECT's surface mesh is a 2D surface that is 3x3 (nx x ny) points,
          // FastScape will be set as a 2D 5x5 point surface. On return to ASPECT, the outer ghost nodes
          // will be ignored, and ASPECT will see only the inner 3x3 surface of FastScape.
          if (use_ghost_nodes)
            set_ghost_nodes(h.get(), vx.get(), vy.get(), vz.get());

          // If specified, apply the orographic controls to the FastScape model.
          if (use_orographic_controls)
            apply_orographic_controls(h.get(), kd.get(), kf.get());

          // Get current FastScape timestep.
          int istep = 0;
          fastscape_get_step_(&istep);

          // Write a file to store h & step in case of restart.
          // TODO: It would be good to roll this isnto the general ASPECT checkpointing,
          // and when we take this needs to be changed.
          if ((this->get_parameters().checkpoint_time_secs == 0) &&
              (this->get_parameters().checkpoint_steps > 0) &&
              (current_timestep % this->get_parameters().checkpoint_steps == 0))
            {

              std::ofstream out_h (restart_filename.c_str());
              std::ofstream out_step (restart_step_filename.c_str());
              std::ofstream out_b (restart_filename_basement.c_str());
              std::stringstream bufferb;
              std::stringstream bufferh;

              fastscape_copy_basement_(b.get());

              out_step << (istep+restart_step) << "\n";

              for (int i=0; i<array_size; i++)
                {
                  bufferh << h[i] << "\n";
                  bufferb << b[i] << "\n";
                }

              out_h << bufferh.str();
              out_b << bufferb.str();
            }

          // Set velocity components.
          if (use_velocities)
            {
              fastscape_set_u_(vz.get());
              fastscape_set_v_(vx.get(), vy.get());
            }

          // Set h to new values, and erosional parameters if there have been changes.
          fastscape_set_h_(h.get());
          fastscape_set_erosional_parameters_(kf.get(), &kfsed, &m, &n, kd.get(), &kdsed, &g, &g, &p);

          // Find  timestep size, run fastscape, and make visualizations.
          execute_fastscape(h.get(), kd.get(), istep);

          // Find out our velocities from the change in height.
          // Where V is a vector of array size that exists on all processes.
          for (int i=0; i<array_size; i++)
            {
              V[i] = (h[i] - h_old[i])/a_dt;
            }

          MPI_Bcast(&V[0], array_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }
      else
        {
          for (unsigned int i=0; i<temporary_variables.size(); i++)
            MPI_Ssend(&temporary_variables[i][0], temporary_variables[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          MPI_Bcast(&V[0], array_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }

      // Get the sizes needed for a data table of the mesh velocities.
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        {
          size_idx[d] = table_intervals[d]+1;
        }

      // Initialize a table to hold all velocity values that will be interpolated back to ASPECT.
      Table<dim,double> velocity_table = fill_data_table(V, size_idx, nx, ny);

      // As our grid_extent variable end points do not account for the change related to an origin
      // not at 0, we adjust this here into an interpolation extent.
      std::array<std::pair<double,double>,dim> interpolation_extent;
      for (unsigned int i=0; i<dim; ++i)
        {
          interpolation_extent[i].first = grid_extent[i].first;
          interpolation_extent[i].second = (grid_extent[i].second + grid_extent[i].first);
        }

      Functions::InterpolatedUniformGridData<dim> *velocities;
      velocities = new Functions::InterpolatedUniformGridData<dim> (interpolation_extent,
                                                                    table_intervals,
                                                                    velocity_table);

      VectorFunctionFromScalarFunctionObject<dim> vector_function_object(
        [&](const Point<dim> &p) -> double
      {
        return velocities->value(p);
      },
      dim-1,
      dim);

      VectorTools::interpolate_boundary_values (mesh_deformation_dof_handler,
                                                *boundary_ids.begin(),
                                                vector_function_object,
                                                mesh_velocity_constraints);
    }


    template <int dim>
    std::vector<std::vector<double>>
    FastScape<dim>::get_aspect_values() const
    {

      const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
      std::vector<std::vector<double>> temporary_variables(dim+2, std::vector<double>());

      // Get a quadrature rule that exists only on the corners, and increase the refinement if specified.
      const QIterated<dim-1> face_corners (QTrapez<1>(),
                                           pow(2,additional_refinement+surface_refinement_difference));

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

                std::vector<Tensor<1,dim>> vel(face_corners.size());
                fe_face_values.reinit(cell, face_no);
                fe_face_values[this->introspection().extractors.velocities].get_function_values(this->get_solution(), vel);

                for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                  {
                    const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                    // Find what x point we're at. Add 1 or 2 depending on if ghost nodes are used.
                    const double indx = 1+use_ghost_nodes+(vertex(0) - grid_extent[0].first)/dx;

                    // If our x or y index isn't close to a whole number, then it's likely an artifact
                    // from using an over-resolved quadrature rule, in that case ignore it.
                    if (abs(indx - round(indx)) >= precision)
                      continue;


                    // If we're in 2D, we want to take the values and apply them to every row of X points.
                    if (dim == 2)
                      {
                        for (int ys=0; ys<ny; ys++)
                          {
                            // FastScape indexes from 1 to n, starting at X and Y = 0, and increases
                            // across the X row. At the end of the row, it jumps back to X = 0
                            // and up to the next X row in increasing Y direction. We track
                            // this to correctly place the variables later on.
                            // Nx*ys effectively tells us what row we are in
                            // and then indx tells us what position in that row.
                            const double index = round(indx)+nx*ys;

                            temporary_variables[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);
                            temporary_variables[1].push_back(index-1);

                            for (unsigned int i=0; i<dim; ++i)
                              {
                                // Always convert to m/yr for FastScape
                                temporary_variables[i+2].push_back(vel[corner][i]*year_in_seconds);
                              }
                          }
                      }
                    // 3D case
                    else
                      {
                        // Because indy only gives us the row we're in, we don't need to add 2 for the ghost node.
                        const double indy = 1+use_ghost_nodes+(vertex(1) - grid_extent[1].first)/dy;

                        if (abs(indy - round(indy)) >= precision)
                          continue;

                        const double index = round((indy-1))*nx+round(indx);

                        temporary_variables[0].push_back(vertex(dim-1) - grid_extent[dim-1].second);   //z component
                        temporary_variables[1].push_back(index-1);

                        for (unsigned int i=0; i<dim; ++i)
                          {
                            temporary_variables[i+2].push_back(vel[corner][i]*year_in_seconds);
                          }
                      }
                  }
              }

      return temporary_variables;
    }


    template <int dim>
    void FastScape<dim>::fill_fastscape_arrays(double *h, double *kd, double *kf, double *vx, double *vy, double *vz, std::vector<std::vector<double>> temporary_variables) const
    {
      // Initialize kf and kd.
      for (int i=0; i<array_size; i++)
        {
          kf[i] = kff;
          kd[i] = kdd;
        }

      for (unsigned int i=0; i<temporary_variables[1].size(); i++)
        {

          int index = temporary_variables[1][i];
          h[index] = temporary_variables[0][i];
          vx[index] = temporary_variables[2][i];
          vz[index] = temporary_variables[dim+1][i];

          if (dim == 2 )
            vy[index] = 0;
          else
            vy[index] = temporary_variables[3][i];
        }

      for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
        {
          // First, find out the size of the array a process wants to send.
          MPI_Status status;
          MPI_Probe(p, 42, this->get_mpi_communicator(), &status);
          int incoming_size = 0;
          MPI_Get_count(&status, MPI_DOUBLE, &incoming_size);

          // Resize the array so it fits whatever the process sends.
          for (unsigned int i=0; i<temporary_variables.size(); ++i)
            {
              temporary_variables[i].resize(incoming_size);
            }

          for (unsigned int i=0; i<temporary_variables.size(); i++)
            MPI_Recv(&temporary_variables[i][0], incoming_size, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

          // Now, place the numbers into the correct place based off the index.
          for (unsigned int i=0; i<temporary_variables[1].size(); i++)
            {
              int index = temporary_variables[1][i];
              h[index] = temporary_variables[0][i];
              vx[index] = temporary_variables[2][i];
              vz[index] = temporary_variables[dim+1][i];

              // In 2D there are no y velocities, so we set them to zero.
              if (dim == 2 )
                vy[index] = 0;
              else
                vy[index] = temporary_variables[3][i];
            }
        }
    }


    template <int dim>
    void FastScape<dim>::initialize_fastscape(double *h, double *b, double *kd, double *kf) const
    {
      const int current_timestep = this->get_timestep_number ();

      // Initialize FastScape with grid and extent.
      fastscape_init_();
      fastscape_set_nx_ny_(&nx,&ny);
      fastscape_setup_();
      fastscape_set_xl_yl_(&x_extent,&y_extent);

      // Set boundary conditions
      fastscape_set_bc_(&bc);

      // Initialize topography
      fastscape_init_h_(h);

      // Set erosional parameters. May have to move this if sed values are updated over time.
      fastscape_set_erosional_parameters_(kf, &kfsed, &m, &n, kd, &kdsed, &g, &gsed, &p);

      if (use_marine)
        fastscape_set_marine_parameters_(&sl, &p1, &p2, &z1, &z2, &r, &l, &kds1, &kds2);

      // Only set the basement if it's a restart
      if (current_timestep != 1)
        fastscape_set_basement_(b);
    }


    template <int dim>
    void FastScape<dim>::execute_fastscape(double *h, double *kd, int istep) const
    {
      const double a_dt = this->get_timestep()/year_in_seconds;

      // Because on the first timestep we will create an initial VTK file before running FastScape
      // and a second after, we first set the visualization step to zero.
      int visualization_step = 0;
      const int current_timestep = this->get_timestep_number ();
      std::string dirname = this->get_output_directory();
      const char *c=dirname.c_str();
      int length = dirname.length();

      // Find a FastScape timestep that is below our maximum timestep.
      int fastscape_iterations = nstep;
      double f_dt = a_dt/fastscape_iterations;
      while (f_dt>maximum_fastscape_timestep)
        {
          fastscape_iterations=fastscape_iterations*2;
          f_dt = a_dt/fastscape_iterations;
        }

      // Set time step
      fastscape_set_dt_(&f_dt);
      fastscape_iterations = fastscape_iterations + istep;
      this->get_pcout() << "   Calling FastScape... " << (fastscape_iterations-istep) << " timesteps of " << f_dt << " years." << std::endl;
      {
        auto t_start = std::chrono::high_resolution_clock::now();

        // If we use stratigraphy it'll handle visualization and not the normal function.
        // TODO: The frequency in this needs to be the same as the total timesteps FastScape will
        // run for, need to figure out how to work this in better.
        if (use_stratigraphy && current_timestep == 1)
          fastscape_strati_(&nstepp, &nreflectorp, &fastscape_iterations, &vexp);
        else if (!use_stratigraphy && current_timestep == 1)
          {
            this->get_pcout() << "      Writing initial VTK..." << std::endl;
            // Note: Here, the HHHHH field in visualization is set to show the diffusivity. However, you can change this so any parameter
            // is visualized.
            fastscape_named_vtk_(kd, &vexp, &visualization_step, c, &length);
          }

        do
          {
            // Execute step, this increases timestep counter
            fastscape_execute_step_();

            // Get value of time step counter
            fastscape_get_step_(&istep);

            // Outputs new h values
            fastscape_copy_h_(h);
          }
        while (istep<fastscape_iterations);

        // Output how long FastScape took to run.
        auto t_end = std::chrono::high_resolution_clock::now();
        double r_time = std::chrono::duration<double>(t_end-t_start).count();
        this->get_pcout() << "      FastScape runtime... " << round(r_time*1000)/1000 << "s" << std::endl;
      }

      visualization_step = current_timestep;

      // Determine whether to create a VTK file this timestep.
      bool make_vtk = 0;
      if (this->get_time() >= last_output_time + output_interval || this->get_time()+this->get_timestep() >= end_time)
        {
          // Don't create a visualization file on a restart.
          if (!restart)
            make_vtk = 1;

          if (output_interval > 0)
            {
              // We need to find the last time output was supposed to be written.
              // this is the last_output_time plus the largest positive multiple
              // of output_intervals that passed since then. We need to handle the
              // edge case where last_output_time+output_interval==current_time,
              // we did an output and std::floor sadly rounds to zero. This is done
              // by forcing std::floor to round 1.0-eps to 1.0.
              const double magic = 1.0+2.0*std::numeric_limits<double>::epsilon();
              last_output_time = last_output_time + std::floor((this->get_time()-last_output_time)/output_interval*magic) * output_interval/magic;
            }
        }

      if (make_vtk)
        {
          this->get_pcout() << "      Writing VTK..." << std::endl;
          fastscape_named_vtk_(kd, &vexp, &visualization_step, c, &length);
        }

      // If we've reached the end time, destroy FastScape.
      if (this->get_time()+this->get_timestep() > end_time)
        {
          this->get_pcout() << "      Destroying FastScape..." << std::endl;
          fastscape_destroy_();
        }
    }


    template <int dim>
    void FastScape<dim>::apply_orographic_controls(double *h, double *kd, double *kf) const
    {
      // First for the wind barrier, we find the maximum height and index
      // along each line in the x and y direction.
      // If wind is east or west, we find maximum point for each ny row along x.
      std::vector<std::vector<double>> hmaxx(2, std::vector<double>(ny, 0.0));
      if (wind_direction == 0 || wind_direction == 1)
        {
          for (int i=0; i<ny; i++)
            {
              for (int j=0; j<nx; j++)
                {
                  if ( h[nx*i+j] > hmaxx[0][i])
                    {
                      hmaxx[0][i] = h[nx*i+j];
                      hmaxx[1][i] = j;
                    }
                }
            }
        }

      // If wind is north or south, we find maximum point for each nx row along y.
      std::vector<std::vector<double>> hmaxy(2, std::vector<double>(nx, 0.0));
      if (wind_direction == 2 || wind_direction == 3)
        {
          for (int i=0; i<nx; i++)
            {
              for (int j=0; j<ny; j++)
                {
                  if ( h[nx*j+i] > hmaxy[0][i])
                    {
                      hmaxy[0][i] = h[nx*j+i];
                      hmaxy[1][i] = j;
                    }
                }
            }
        }

      // Now we loop through all the points again and apply the factors.
      // TODO: I made quite a few changes here, should double check it still works in all cases.
      std::vector<double> control_applied(array_size, 0);
      for (int i=0; i<ny; i++)
        {
          // Factor from wind barrier. Apply a switch based off wind direction.
          // Where 0 is wind going to the west, 1 the east, 2 the south, and 3 the north.
          for (int j=0; j<nx; j++)
            {
              switch (wind_direction)
                {
                  case 0 :
                  {
                    // If we are above the set elevation, and on the correct side based on the wind direction apply
                    // the factor. Apply this regardless of whether or not we stack controls.
                    if ( (hmaxx[0][i] > wind_barrier_elevation) && (j < hmaxx[1][i]) )
                      {
                        kf[nx*i+j] = kf[nx*i+j]*wind_barrier_erosional_factor;
                        kd[nx*i+j] = kd[nx*i+j]*wind_barrier_erosional_factor;
                        control_applied[nx*i+j] = 1;
                      }
                    break;
                  }
                  case 1 :
                  {
                    if ( (hmaxx[0][i] > wind_barrier_elevation) && (j > hmaxx[1][i]) )
                      {
                        kf[nx*i+j] = kf[nx*i+j]*wind_barrier_erosional_factor;
                        kd[nx*i+j] = kd[nx*i+j]*wind_barrier_erosional_factor;
                        control_applied[nx*i+j] = 1;
                      }
                    break;
                  }
                  case 2 :
                  {
                    if ( (hmaxy[0][j] > wind_barrier_elevation) && (i > hmaxy[1][j]) )
                      {
                        kf[nx*i+j] = kf[nx*i+j]*wind_barrier_erosional_factor;
                        kd[nx*i+j] = kd[nx*i+j]*wind_barrier_erosional_factor;
                        control_applied[nx*i+j] = 1;
                      }
                    break;
                  }
                  case 3 :
                  {
                    if ( (hmaxy[0][j] > wind_barrier_elevation) && (i < hmaxy[1][j]) )
                      {
                        kf[nx*i+j] = kf[nx*i+j]*wind_barrier_erosional_factor;
                        kd[nx*i+j] = kd[nx*i+j]*wind_barrier_erosional_factor;
                        control_applied[nx*i+j] = 1;
                      }
                    break;
                  }
                  default :
                    AssertThrow(false, ExcMessage("This does not correspond with a wind direction."));
                    break;
                }

              // If we are above the flat elevation and stack controls, apply the flat elevation factor. If we are not
              // stacking controls, apply the factor if the wind barrier was not applied to this point.
              if ( ((h[nx*i+j] > flat_elevation) && stack_controls==true) || ((h[nx*i+j] > flat_elevation) && !stack_controls && (control_applied[nx*i+j]==0)) )
                {
                  kf[nx*i+j] = kf[nx*i+j]*flat_erosional_factor;
                  kd[nx*i+j] = kd[nx*i+j]*flat_erosional_factor;
                }
              // If we are not stacking controls and the wind barrier was applied to this point, only
              // switch to this control if the factor is greater.
              else if ( (h[nx*i+j] > flat_elevation) && stack_controls==false && (control_applied[nx*i+j]==1) && (flat_erosional_factor > wind_barrier_erosional_factor) )
                {
                  if ( wind_barrier_erosional_factor != 0)
                    {
                      kf[nx*i+j] = (kf[nx*i+j]/wind_barrier_erosional_factor)*flat_erosional_factor;
                      kd[nx*i+j] = (kd[nx*i+j]/wind_barrier_erosional_factor)*flat_erosional_factor;
                    }
                  // If a wind barrier factor of zero was applied for some reason, we set it back to the default
                  // and apply the flat_erosional_factor.
                  else
                    {
                      kf[nx*i+j] = kff*flat_erosional_factor;
                      kd[nx*i+j] = kdd*flat_erosional_factor;
                    }
                }
            }
        }
    }


    template <int dim>
    void FastScape<dim>::set_ghost_nodes(double *h, double *vx, double *vy, double *vz) const
    {
      const int current_timestep = this->get_timestep_number ();
      std::unique_ptr<double[]> slopep (new double[array_size]());

      // Copy the slopes at each point, this will be used to set an H
      // at the ghost nodes if a boundary mass flux is given.
      fastscape_copy_slope_(slopep.get());

      // Here we set the ghost nodes at the left and right boundaries. In most cases,
      // this involves setting the node to the same values of v and h as the inward node.
      // With the inward node being above or below for the bottom and top rows of ghost nodes,
      // or to the left and right for the right and left columns of ghost nodes.
      for (int j=0; j<ny; j++)
        {
          // Nx*j will give us the row we're in, and one is subtracted as FastScape starts from 1 not zero.
          // If we're on the left, the multiple of the row will always represent the first node.
          // Subtracting one from the row above this gives us the last node of the previous row.
          const int index_left = nx*j;
          const int index_right = nx*(j+1)-1;
          double slope = 0;

          // Here we set the ghost nodes to the value of the nodes next to them, where for the left we
          // add one to go to the node to the right, and for the right side
          // we subtract one to go to the inner node to the left.
          vz[index_right] = vz[index_right-1];
          vz[index_left] =  vz[index_left+1];

          vy[index_right] = vy[index_right-1];
          vy[index_left] = vy[index_left+1];

          vx[index_right] = vx[index_right-1];
          vx[index_left] = vx[index_left+1];

          if (current_timestep == 1 || left_flux == 0)
            {
              // If it's the first timestep add in initial slope. If we have no flux,
              // set the ghost node to the node next to it.
              // FastScape calculates the slope by looking at all nodes surrounding the point
              // so we need to consider the slope over 2 dx.
              slope = left_flux/kdd;
              h[index_left] = h[index_left+1] + slope*2*dx;

              // This sets the ghost nodes of fixed boundaries to a user specified height which will
              // be used as an erosional baselevel by fastscape and may limit erosional flux out
              // of the model domain
              if (use_extra_base_level)
                {
                  if ((left==1) && (h[index_left]<h_extra_base_level))
                    h[index_left] = h_extra_base_level;
                }
            }
          else
            {
              // If we have flux through a boundary, we need to update the height to keep the correct slope.
              // Because the corner nodes always show a slope of zero, this will update them according to
              // dthe closest non-ghost node. E.g. if we're at a corner node, look instead up a row and inward.
              if (j == 0)
                slope = left_flux/kdd - std::tan(slopep[index_left+nx+1]*numbers::PI/180.);
              else if (j==(ny-1))
                slope = left_flux/kdd - std::tan(slopep[index_left-nx+1]*numbers::PI/180.);
              else
                slope = left_flux/kdd - std::tan(slopep[index_left+1]*numbers::PI/180.);

              h[index_left] = h[index_left] + slope*2*dx;
            }

          if (current_timestep == 1 || right_flux == 0)
            {
              slope = right_flux/kdd;
              h[index_right] = h[index_right-1] + slope*2*dx;
              if (use_extra_base_level)
                {
                  if ((right==1) && (h[index_right]<h_extra_base_level))
                    h[index_right] = h_extra_base_level;
                }
            }
          else
            {
              if (j == 0)
                slope = right_flux/kdd - std::tan(slopep[index_right+nx-1]*numbers::PI/180.);
              else if (j==(ny-1))
                slope = right_flux/kdd - std::tan(slopep[index_right-nx-1]*numbers::PI/180.);
              else
                slope = right_flux/kdd - std::tan(slopep[index_right-1]*numbers::PI/180.);

              h[index_right] = h[index_right] + slope*2*dx;
            }

          // This sets the ghost nodes of fixed boundaries to a user specified height which will
          // be used as an erosional baselevel by fastscape and may limit erosional flux out
          // of the model domain
          if (use_extra_base_level)
            {
              if (right==1
                  && h[index_right]<h_extra_base_level
                  && leftright_ghost_nodes_periodic == false)
                h[index_right] = h_extra_base_level;
              if (left==1
                  && h[index_left]<h_extra_base_level
                  && leftright_ghost_nodes_periodic == false)
                h[index_left] = h_extra_base_level;
            }

          // If the boundaries are periodic, then we look at the velocities on both sides of the
          // model, and set the ghost node according to the direction of flow. As FastScape will
          // receive all velocities it will have a direction, and we only need to look at the (non-ghost)
          // nodes directly to the left and right.
          if (left == 0 && right == 0 || leftright_ghost_nodes_periodic == true)
            {
              // First we assume that flow is going to the left.
              int side = index_left;
              int op_side = index_right;

              // Indexing depending on which side the ghost node is being set to.
              int jj = 1;

              // If nodes on both sides are going the same direction, then set the respective
              // ghost nodes to equal these sides. By doing this, the ghost nodes at the opposite
              // side of flow will work as a mirror mimicking what is happening on the other side.
              if (vx[index_right-1] > 0 && vx[index_left+1] >= 0)
                {
                  side = index_right;
                  op_side = index_left;
                  jj = -1;
                }
              else if (vx[index_right-1] <= 0 && vx[index_left+1] < 0)
                {
                  side = index_left;
                  op_side = index_right;
                  jj = 1;
                }
              else
                continue;

              // Set right ghost node
              h[index_right] = h[side+jj];
              vx[index_right] = vx[side+jj];
              vy[index_right] = vy[side+jj];
              vz[index_right] = vz[side+jj];

              // Set left ghost node
              h[index_left] = h[side+jj];
              vx[index_left] = vx[side+jj];
              vy[index_left] = vy[side+jj];
              vz[index_left] = vz[side+jj];

              // Set opposing ASPECT boundary so it's periodic.
              h[op_side-jj] = h[side+jj];
              vx[op_side-jj] = vx[side+jj];
              vz[op_side-jj] = vz[side+jj];
              vy[op_side-jj] = vy[side+jj];
            }
        }

      // Now do the same for the top and bottom ghost nodes.
      for (int j=0; j<nx; j++)
        {
          // The bottom row indexes are 0 to nx-1.
          const int index_bot = j;

          // Nx multiplied by (total rows - 1) gives us the start of
          // the top row, and j gives the position in the row.
          const int index_top = nx*(ny-1)+j;
          double slope = 0;

          vz[index_bot] = vz[index_bot+nx];
          vz[index_top] = vz[index_top-nx];

          vy[index_bot] = vy[index_bot+nx];
          vy[index_top] = vy[index_top-nx];

          vx[index_bot] = vx[index_bot+nx];
          vx[index_top] = vx[index_top-nx];

          if (current_timestep == 1 || top_flux == 0)
            {
              slope = top_flux/kdd;
              h[index_top] = h[index_top-nx] + slope*2*dx;
              if (use_extra_base_level)
                {
                  if ((top==1) && (h[index_top]<h_extra_base_level))
                    h[index_top] = h_extra_base_level;
                }
            }
          else
            {
              if (j == 0)
                slope = top_flux/kdd - std::tan(slopep[index_top-nx+1]*numbers::PI/180.);
              else if (j==(nx-1))
                slope = top_flux/kdd - std::tan(slopep[index_top-nx-1]*numbers::PI/180.);
              else
                slope = top_flux/kdd - std::tan(slopep[index_top-nx]*numbers::PI/180.);

              h[index_top] = h[index_top] + slope*2*dx;
              if (use_extra_base_level)
                {
                  if ((bottom==1) && (h[index_bot]<h_extra_base_level))
                    h[index_bot] = h_extra_base_level;
                }
            }

          if (current_timestep == 1 || bottom_flux == 0)
            {
              slope = bottom_flux/kdd;
              h[index_bot] = h[index_bot+nx] + slope*2*dx;
            }
          else
            {
              if (j == 0)
                slope = bottom_flux/kdd - std::tan(slopep[index_bot+nx+1]*numbers::PI/180.);
              else if (j==(nx-1))
                slope = bottom_flux/kdd - std::tan(slopep[index_bot+nx-1]*numbers::PI/180.);
              else
                slope = bottom_flux/kdd - std::tan(slopep[index_bot+nx]*numbers::PI/180.);

              h[index_bot] = h[index_bot] + slope*2*dx;
            }

          if (use_extra_base_level)
            {
              if (bottom==1
                  && h[index_bot]<h_extra_base_level
                  && topbottom_ghost_nodes_periodic == false)
                h[index_bot] = h_extra_base_level;
              if (top==1
                  && h[index_top]<h_extra_base_level
                  && topbottom_ghost_nodes_periodic == false)
                h[index_top] = h_extra_base_level;
            }

          if (bottom == 0 && top == 0 || topbottom_ghost_nodes_periodic == true)
            {
              int side = index_bot;
              int op_side = index_top;
              int jj = nx;

              if (vy[index_bot+nx-1] > 0 && vy[index_top-nx-1] >= 0)
                {
                  side = index_top;
                  op_side = index_bot;
                  jj = -nx;
                }
              else if (vy[index_bot+nx-1] <= 0 && vy[index_top-nx-1] < 0)
                {
                  side = index_bot;
                  op_side = index_top;
                  jj = nx;
                }
              else
                continue;

              // Set top ghost node
              h[index_top] = h[side+jj];
              vx[index_top] = vx[side+jj];
              vy[index_top] = vy[side+jj];
              vz[index_top] = vz[side+jj];

              // Set bottom ghost node
              h[index_bot] = h[side+jj];
              vx[index_bot] = vx[side+jj];
              vy[index_bot] = vy[side+jj];
              vz[index_bot] = vz[side+jj];

              // Set opposing ASPECT boundary so it's periodic.
              h[op_side-jj] = h[side+jj];
              vx[op_side-jj] = vx[side+jj];
              vz[op_side-jj] = vz[side+jj];
              vy[op_side-jj] = vy[side+jj];
            }
        }
    }


    template <int dim>
    Table<dim,double>
    FastScape<dim>::fill_data_table(std::vector<double> values, TableIndices<dim> size_idx, int nx, int ny) const
    {
      // Create data table based off of the given size.
      Table<dim,double> data_table;
      data_table.TableBase<dim,double>::reinit(size_idx);
      TableIndices<dim> idx;

      // Loop through the data table and fill it with the velocities from FastScape.
      if (dim == 2)
        {
          std::vector<double> V2(nx);

          for (int i=1; i<(nx-1); i++)
            {
              // If using the center slice, find velocities from the row closest to the center.
              if (center_slice)
                {
                  const int index = i+nx*(round((ny-1)/2));
                  V2[i-1] = values[index];
                }
              // Here we use average velocities across the y nodes, excluding the ghost nodes (top and bottom row).
              // Note: If ghost nodes are turned off, boundary effects may influence this.
              else
                {
                  for (int ys=1; ys<(ny-1); ys++)
                    {
                      const int index = i+nx*ys;
                      V2[i-1] += values[index];
                    }
                  V2[i-1] = V2[i-1]/(ny-2);
                }
            }

          for (unsigned int i=0; i<data_table.size()[1]; ++i)
            {
              idx[1] = i;

              for (unsigned int j=0; j<(data_table.size()[0]); ++j)
                {
                  idx[0] = j;

                  // Convert back to m/s.
                  data_table(idx) = V2[j]/year_in_seconds;
                }
            }
        }
      else
        {
          // Indexes through z, y, and then x.
          for (unsigned int k=0; k<data_table.size()[2]; ++k)
            {
              idx[2] = k;

              for (unsigned int i=0; i<data_table.size()[1]; ++i)
                {
                  idx[1] = i;

                  for (unsigned int j=0; j<data_table.size()[0]; ++j)
                    {
                      idx[0] = j;

                      // Convert back to m/s.
                      data_table(idx) = values[(nx+1)*use_ghost_nodes+nx*i+j]/year_in_seconds;

                    }
                }
            }
        }

      return data_table;
    }


    template <int dim>
    void FastScape<dim>::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection ("Fastscape");
        {
          prm.declare_entry("Number of steps", "10",
                            Patterns::Integer(),
                            "Number of steps per ASPECT timestep");
          prm.declare_entry("Maximum timestep", "10e3",
                            Patterns::Double(0),
                            "Maximum timestep for FastScape. Units: $\\{yrs}$");
          prm.declare_entry("Vertical exaggeration", "-1",
                            Patterns::Double(),
                            "Vertical exaggeration for FastScape's VTK file. -1 outputs topography, basement, and sealevel.");
          prm.declare_entry("Additional fastscape refinement", "0",
                            Patterns::Integer(),
                            "How many levels above ASPECT FastScape should be refined.");
          prm.declare_entry ("Use center slice for 2d", "false",
                             Patterns::Bool (),
                             "If this is set to true, then a 2D model will only consider the "
                             "center slice FastScape gives. If set to false, then aspect will"
                             "average the inner third of what FastScape calculates.");
          prm.declare_entry("Fastscape seed", "1000",
                            Patterns::Integer(),
                            "Seed used for adding an initial noise to FastScape topography based on the initial noise magnitude.");
          prm.declare_entry("Maximum surface refinement level", "1",
                            Patterns::Integer(),
                            "This should be set to the highest ASPECT refinement level expected at the surface.");
          prm.declare_entry("Surface refinement difference", "0",
                            Patterns::Integer(),
                            "The difference between the lowest and highest refinement level at the surface. E.g., if three resolution "
                            "levels are expected, this would be set to 2.");
          prm.declare_entry ("Use marine component", "false",
                             Patterns::Bool (),
                             "Flag to use the marine component of FastScape.");
          prm.declare_entry ("Use stratigraphy", "false",
                             Patterns::Bool (),
                             "Flag to use stratigraphy");
          prm.declare_entry("Total steps", "100000",
                            Patterns::Integer(),
                            "Total number of steps you expect in the FastScape model, only used for stratigraphy.");
          prm.declare_entry("Number of horizons", "1",
                            Patterns::Integer(),
                            "Number of horizons to track and visualize in FastScape, only used for stratigraphy");
          prm.declare_entry("Y extent in 2d", "100000",
                            Patterns::Double(),
                            "FastScape Y extent when using a 2D ASPECT model. Units: $\\{m}$");
          prm.declare_entry ("Use ghost nodes", "true",
                             Patterns::Bool (),
                             "Flag to use ghost nodes");
          prm.declare_entry ("Use velocities", "true",
                             Patterns::Bool (),
                             "Flag to use FastScape advection and uplift.");
          prm.declare_entry("Precision", "0.001",
                            Patterns::Double(),
                            "Precision value for how close a ASPECT node must be to the FastScape node for the value to be transferred.");
          prm.declare_entry ("Sediment rain rates", "0,0",
                             Patterns::List (Patterns::Double(0)),
                             "Sediment rain rates given as a list equal to the number of intervals. Units: $\\{m/yr}$");
          prm.declare_entry ("Sediment rain time intervals", "0",
                             Patterns::List (Patterns::Double(0)),
                             "A list of sediment rain times. Units: $\\{yrs}$");
          prm.declare_entry("Initial noise magnitude", "5",
                            Patterns::Integer(),
                            "Maximum topography change from the initial noise. Units: $\\{m}$");

          prm.enter_subsection ("Boundary conditions");
          {
            prm.declare_entry ("Bottom", "1",
                               Patterns::Integer (0, 1),
                               "Bottom boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Right", "1",
                               Patterns::Integer (0, 1),
                               "Right boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Top", "1",
                               Patterns::Integer (0, 1),
                               "Top boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry ("Left", "1",
                               Patterns::Integer (0, 1),
                               "Left boundary condition, where 1 is fixed and 0 is reflective.");
            prm.declare_entry("Left mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through left boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Right mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through right boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Top mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through top boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Bottom mass flux", "0",
                              Patterns::Double(),
                              "Flux per unit length through bottom boundary. Units: $\\{m^2/yr}$ ");
            prm.declare_entry ("Top bottom ghost nodes periodic", "false",
                               Patterns::Bool (),
                               "Whether to set the ghost nodes at the FastScape top and bottom boundary "
                               "to periodic even if 'Bottom' and 'Top' are set to fixed boundary.");
            prm.declare_entry ("Left right ghost nodes periodic", "false",
                               Patterns::Bool (),
                               "Whether to set the ghost nodes at the FastScape left and right boundary "
                               "to periodic even if 'Left' and 'Right' are set to fixed boundary.");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Erosional parameters");
          {
            prm.declare_entry("Drainage area exponent", "0.4",
                              Patterns::Double(),
                              "Exponent for drainage area.");
            prm.declare_entry("Slope exponent", "1",
                              Patterns::Double(),
                              "The  slope  exponent  for  SPL (n).  Generally  m/n  should  equal  approximately 0.4");
            prm.declare_entry("Multi-direction slope exponent", "1",
                              Patterns::Double(),
                              "Exponent to determine the distribution from the SPL to neighbor nodes, with"
                              "10 being steepest decent and 1 being more varied.");
            prm.declare_entry("Bedrock deposition coefficient", "-1",
                              Patterns::Double(),
                              "Deposition coefficient for bedrock");
            prm.declare_entry("Sediment deposition coefficient", "-1",
                              Patterns::Double(),
                              "Deposition coefficient for sediment");
            prm.declare_entry("Bedrock river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for bedrock in the Stream Power Law. Units: $\\{m^(1-2*drainage_area_exponent)/yr}$");
            prm.declare_entry("Sediment river incision rate", "-1",
                              Patterns::Double(),
                              "River incision rate for sediment in the Stream Power Law. Units: $\\{m^(1-2*drainage_area_exponent)/yr}$ ");
            prm.declare_entry("Bedrock diffusivity", "1",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for bedrock. Units: $\\{m^2/yr}$ ");
            prm.declare_entry("Sediment diffusivity", "-1",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for sediment. Units: $\\{m^2/yr}$");
            prm.declare_entry("Orographic elevation control", "2000",
                              Patterns::Integer(),
                              "Above this height, the elevation factor is applied. Units: $\\{m}$");
            prm.declare_entry("Orographic wind barrier height", "500",
                              Patterns::Integer(),
                              "When terrain reaches this height the wind barrier factor is applied. Units: $\\{m}$");
            prm.declare_entry("Elevation factor", "1",
                              Patterns::Double(),
                              "Amount to multiply kf and kd by past given orographic elevation control.");
            prm.declare_entry("Wind barrier factor", "1",
                              Patterns::Double(),
                              "Amount to multiply kf and kd by past given wind barrier height.");
            prm.declare_entry ("Stack orographic controls", "true",
                               Patterns::Bool (),
                               "Whether or not to apply both controls to a point, or only a maximum of one set as the wind barrier.");
            prm.declare_entry ("Flag to use orographic controls", "false",
                               Patterns::Bool (),
                               "Whether or not to apply orographic controls.");
            prm.declare_entry ("Wind direction", "west",
                               Patterns::Selection("east|west|south|north"),
                               "This parameter assumes a wind direction, deciding which side is reduced from the wind barrier.");
            prm.declare_entry ("Use an erosional base level differing from sea level", "false",
                               Patterns::Bool (),
                               "Whether or not to use an erosional base level that differs from sea level. Setting this parameter to "
                               "true will set all ghost nodes of fixed FastScape boundaries to the height you specify in "
                               "'set Erosional base level differing from sea level'. \nThis can make "
                               "sense for a continental model where the model surrounding topography is assumed above sea level, "
                               "e.g. highlands. If the sea level would be used as an erosional base level in this case, all topography "
                               "erodes away with lots of 'sediment volume' lost through the sides of the model. This is mostly "
                               "important, when there are mountains in the middle of the model, while it is less important when there "
                               "is lower relief in the middle of the model. \n"
                               "In the FastScape  visualization files, setting the extra base level may show up as a strong "
                               "slope at the fixed boundaries of the model. However, in the ASPECT visualization files it will not "
                               "show up, as the ghost nodes only exist in FastScape.");
            prm.declare_entry("Erosional base level differing from sea level", "0",
                              Patterns::Double(),
                              "When 'Use an erosional base level differing from sea level' is set to true, all ghost nodes of fixed "
                              "FastScape boundaries where no mass flux is specified by the user (FastScape boundary condition set to 1 "
                              "and 'Left/Right/Bottom/Top mass flux' set to 0) will be fixed to this elevation. The "
                              "reflecting boundaries (FastScape boundary condition set to 0) will not be affected, nor are the "
                              "boundaries where a mass flux is specified. \n"
                              "Units: m");
          }
          prm.leave_subsection();

          prm.enter_subsection ("Marine parameters");
          {
            prm.declare_entry("Sea level", "0",
                              Patterns::Double(),
                              "Sea level relative to the ASPECT surface, where the maximum Z or Y extent in ASPECT is a sea level of zero. Units: $\\{m}$ ");
            prm.declare_entry("Sand porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of sand. ");
            prm.declare_entry("Shale porosity", "0.0",
                              Patterns::Double(),
                              "Porosity of shale. ");
            prm.declare_entry("Sand e-folding depth", "1e3",
                              Patterns::Double(),
                              "E-folding depth for the exponential of the sand porosity law. Units: $\\{m}$");
            prm.declare_entry("Shale e-folding depth", "1e3",
                              Patterns::Double(),
                              "E-folding depth for the exponential of the shale porosity law. Units: $\\{m}$");
            prm.declare_entry("Sand-shale ratio", "0.5",
                              Patterns::Double(),
                              "Ratio of sand to shale for material leaving continent.");
            prm.declare_entry("Depth averaging thickness", "1e2",
                              Patterns::Double(),
                              "Depth averaging for the sand-shale equation. Units: $\\{m}$");
            prm.declare_entry("Sand transport coefficient", "5e2",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for sand. Units: $\\{m^2/yr}$");
            prm.declare_entry("Shale transport coefficient", "2.5e2",
                              Patterns::Double(),
                              "Transport coefficient (diffusivity) for shale. Units: $\\{m^2/yr}$ ");
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
      end_time = prm.get_double ("End time");
      if (prm.get_bool ("Use years in output instead of seconds") == true)
        end_time *= year_in_seconds;
      prm.enter_subsection ("Mesh deformation");
      {
        prm.enter_subsection("Fastscape");
        {
          nstep = prm.get_integer("Number of steps");
          maximum_fastscape_timestep = prm.get_double("Maximum timestep");
          vexp = prm.get_double("Vertical exaggeration");
          additional_refinement = prm.get_integer("Additional fastscape refinement");
          center_slice = prm.get_bool("Use center slice for 2d");
          fs_seed = prm.get_integer("Fastscape seed");
          maximum_surface_refinement_level = prm.get_integer("Maximum surface refinement level");
          surface_refinement_difference = prm.get_integer("Surface refinement difference");
          use_marine = prm.get_bool("Use marine component");
          use_stratigraphy = prm.get_bool("Use stratigraphy");
          nstepp = prm.get_integer("Total steps");
          nreflectorp = prm.get_integer("Number of horizons");
          y_extent_2d = prm.get_double("Y extent in 2d");
          use_ghost_nodes = prm.get_bool("Use ghost nodes");
          use_velocities = prm.get_bool("Use velocities");
          precision = prm.get_double("Precision");
          noise_h = prm.get_integer("Initial noise magnitude");
          sediment_rain_rates = Utilities::string_to_double
                                (Utilities::split_string_list(prm.get ("Sediment rain rates")));
          sediment_rain_times = Utilities::string_to_double
                                (Utilities::split_string_list(prm.get ("Sediment rain time intervals")));

          if (!this->convert_output_to_years())
            {
              maximum_fastscape_timestep /= year_in_seconds;
              for (unsigned int j=0; j<sediment_rain_rates.size(); j++)
                sediment_rain_rates[j] *= year_in_seconds;
            }

          if (sediment_rain_rates.size() != sediment_rain_times.size()+1)
            AssertThrow(false, ExcMessage("Error: There must be one more sediment rain value than interval."));

          prm.enter_subsection("Boundary conditions");
          {
            bottom = prm.get_integer("Bottom");
            right = prm.get_integer("Right");
            top = prm.get_integer("Top");
            left = prm.get_integer("Left");
            left_flux = prm.get_double("Left mass flux");
            right_flux = prm.get_double("Right mass flux");
            top_flux = prm.get_double("Top mass flux");
            bottom_flux = prm.get_double("Bottom mass flux");

            if (!this->convert_output_to_years())
              {
                left_flux *= year_in_seconds;
                right_flux *= year_in_seconds;
                top_flux *= year_in_seconds;
                bottom_flux *= year_in_seconds;
              }

            // Put the boundary condition values into a four digit value to send to FastScape.
            bc = bottom*1000+right*100+top*10+left;

            if ((left_flux != 0 && top_flux != 0) || (left_flux != 0 && bottom_flux != 0) ||
                (right_flux != 0 && bottom_flux != 0) || (right_flux != 0 && top_flux != 0))
              AssertThrow(false,ExcMessage("Currently the plugin does not support mass flux through adjacent boundaries."));

            topbottom_ghost_nodes_periodic = prm.get_bool("Top bottom ghost nodes periodic");
            leftright_ghost_nodes_periodic = prm.get_bool("Left right ghost nodes periodic");
          }
          prm.leave_subsection();

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
            flat_elevation = prm.get_integer("Orographic elevation control");
            wind_barrier_elevation = prm.get_integer("Orographic wind barrier height");
            flat_erosional_factor = prm.get_double("Elevation factor");
            wind_barrier_erosional_factor = prm.get_double("Wind barrier factor");
            stack_controls = prm.get_bool("Stack orographic controls");
            use_orographic_controls = prm.get_bool("Flag to use orographic controls");

            if (!this->convert_output_to_years())
              {
                kff *= year_in_seconds;
                kdd *= year_in_seconds;
                kfsed *= year_in_seconds;
                kdsed *= year_in_seconds;
              }

            // Wind direction
            if (prm.get ("Wind direction") == "west")
              wind_direction = 0;
            else if (prm.get ("Wind direction") == "east")
              wind_direction = 1;
            else if (prm.get ("Wind direction") == "north")
              wind_direction = 2;
            else if (prm.get ("Wind direction") == "south")
              wind_direction = 3;
            else
              AssertThrow(false, ExcMessage("Not a valid wind direction."));

            // set fixed ghost nodes to a base level for erosion that differs from sea level
            use_extra_base_level = prm.get_bool("Use an erosional base level differing from sea level");
            if (use_extra_base_level)
              AssertThrow(use_extra_base_level && use_ghost, ExcMessage(
                            "If you want to use an erosional base level differing from sea level, "
                            "you need to use ghost nodes."));
            h_extra_base_level = prm.get_double("Erosional base level differing from sea level");
          }
          prm.leave_subsection();

          prm.enter_subsection("Marine parameters");
          {
            sl = prm.get_double("Sea level");
            p1 = prm.get_double("Sand porosity");
            p2 = prm.get_double("Shale porosity");
            z1 = prm.get_double("Sand e-folding depth");
            z2 = prm.get_double("Shale e-folding depth");
            r = prm.get_double("Sand-shale ratio");
            l = prm.get_double("Depth averaging thickness");
            kds1 = prm.get_double("Sand transport coefficient");
            kds2 = prm.get_double("Shale transport coefficient");

            if (!this->convert_output_to_years())
              {
                kds1 *= year_in_seconds;
                kds2 *= year_in_seconds;
              }
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
      prm.leave_subsection ();

      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Visualization");
        {
          output_interval = prm.get_double ("Time between graphical output");
          if (this->convert_output_to_years())
            output_interval *= year_in_seconds;
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
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
                                           "A plugin that uses the program FastScape to add surface processes "
                                           "such as diffusion, sedimentation, and the Stream Power Law to ASPECT. "
                                           "FastScape is initialized with ASPECT height and velocities, and then "
                                           "continues to run in the background, updating with new ASPECT velocities "
                                           "when called. Once FastScape has run for the given amount of timesteps per ASPECT timestep."
                                           "This plugin compares the initial and final heights to compute a vertical mesh deformation "
                                           "velocity at the surface boundary. More information on FastScape can be found at: https://fastscape.org/fastscapelib-fortran/")
  }
}
#endif
