/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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


#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/initial_topography_model/zero_topography.h>
#include <aspect/simulator_signals.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_tools.h>

#include <numeric>


namespace aspect
{
  namespace GeometryModel
  {
    template <int dim>
    void
    Box<dim>::initialize ()
    {
      // Get pointer to initial topography model
      topo_model = const_cast<InitialTopographyModel::Interface<dim>*>(&this->get_initial_topography_model());

      //surface_height.resize(dim);

      // Check that initial topography is required.
      // If so, connect the initial topography function
      // to the right signals: It should be applied after
      // the final initial adaptive refinement and after a restart.
      if (Plugins::plugin_type_matches<InitialTopographyModel::ZeroTopography<dim>>(*topo_model) == false)
        {
          this->get_signals().pre_set_initial_state.connect(
            [&](typename parallel::distributed::Triangulation<dim> &tria)
          {
            this->topography(tria);
          }
          );
          this->get_signals().post_resume_load_user_data.connect(
            [&](typename parallel::distributed::Triangulation<dim> &tria)
          {
            this->topography(tria);
          }
          );
        }
    }


    template <int dim>
    void
    Box<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      const std::vector<unsigned int> rep_vec(repetitions.begin(), repetitions.end());
      GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                                 rep_vec,
                                                 box_origin,
                                                 box_origin+extents,
                                                 true);

      // Tell p4est about the periodicity of the mesh.
      std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator>>
      periodicity_vector;
      for (int i=0; i<dim; ++i)
        if (periodic[i])
          GridTools::collect_periodic_faces
          ( coarse_grid, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
            /*direction*/ i, periodicity_vector);

      if (periodicity_vector.size() > 0)
        coarse_grid.add_periodicity (periodicity_vector);
    }

    template <int dim>
    void
    Box<dim>::
    topography (typename parallel::distributed::Triangulation<dim> &grid) const
    {
      // Here we provide GridTools with the function to displace vertices
      // in the vertical direction by an amount specified by the initial topography model
      GridTools::transform(
        [&](const Point<dim> &p) -> Point<dim>
      {
        return this->add_topography(p);
      },
      grid);

      this->get_pcout() << "   Added initial topography to grid" << std::endl << std::endl;
    }


    template <int dim>
    Point<dim>
    Box<dim>::
    add_topography (const Point<dim> &x_y_z) const
    {
      // Get the surface x (,y) point
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = x_y_z[d];

      // Get the surface topography at this point
      const double topo = topo_model->value(surface_point);

      // Compute the displacement of the z coordinate
      const double ztopo = (x_y_z[dim-1] - box_origin[dim-1]) / extents[dim-1] * topo;

      // Compute the new point
      Point<dim> x_y_ztopo = x_y_z;
      x_y_ztopo[dim-1] += ztopo;

      return x_y_ztopo;
    }


    template <int dim>
    std::set<types::boundary_id>
    Box<dim>::
    get_used_boundary_indicators () const
    {
      // boundary indicators are zero through 2*dim-1
      std::set<types::boundary_id> s;
      for (unsigned int i=0; i<2*dim; ++i)
        s.insert (i);
      return s;
    }


    template <int dim>
    std::map<std::string,types::boundary_id>
    Box<dim>::
    get_symbolic_boundary_names_map () const
    {
      switch (dim)
        {
          case 2:
          {
            return
            {
              {"left",   0},
              {"right",  1},
              {"bottom", 2},
              {"top",    3}
            };
          }

          case 3:
          {
            return
            {
              {"left",   0},
              {"right",  1},
              {"front",  2},
              {"back",   3},
              {"bottom", 4},
              {"top",    5}

            };
          }
        }

      Assert (false, ExcNotImplemented());
      return {};
    }


    template <int dim>
    std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>>
    Box<dim>::
    get_periodic_boundary_pairs () const
    {
      std::set<std::pair<std::pair<types::boundary_id, types::boundary_id>, unsigned int>> periodic_boundaries;
      for ( unsigned int i=0; i<dim; ++i)
        if (periodic[i])
          periodic_boundaries.insert( std::make_pair( std::pair<types::boundary_id, types::boundary_id>(2*i, 2*i+1), i) );
      return periodic_boundaries;
    }



    template <int dim>
    void
    Box<dim>::adjust_positions_for_periodicity (Point<dim> &position,
                                                const ArrayView<Point<dim>> &connected_positions,
                                                const ArrayView<Tensor<1, dim>> &/*connected_velocities*/) const
    {
      for (unsigned int i = 0; i < dim; ++i)
        if (periodic[i])
          {
            if (position[i] < box_origin[i])
              {
                position[i] += extents[i];
                for (auto &connected_position: connected_positions)
                  connected_position[i] += extents[i];
              }
            else if (position[i] > box_origin[i] + extents[i])
              {
                position[i] -= extents[i];
                for (auto &connected_position: connected_positions)
                  connected_position[i] -= extents[i];
              }
          }
    }



    template <int dim>
    Point<dim>
    Box<dim>::get_extents () const
    {
      return extents;
    }

    template <int dim>
    const std::array<unsigned int, dim> &
    Box<dim>::get_repetitions () const
    {
      return repetitions;
    }

    template <int dim>
    Point<dim>
    Box<dim>::get_origin () const
    {
      return box_origin;
    }

    template <int dim>
    double
    Box<dim>::
    length_scale () const
    {
      return 0.01*extents[0];
    }


    template <int dim>
    double
    Box<dim>::depth(const Point<dim> &position) const
    {
      // Get the surface x (,y) point
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = position[d];

      // Get the surface topography at this point
      const double topo = topo_model->value(surface_point);

      const double d = extents[dim-1] + topo - (position(dim-1)-box_origin[dim-1]);

      return std::min (std::max (d, 0.), maximal_depth());

    }


    template <int dim>
    double
    Box<dim>::height_above_reference_surface(const Point<dim> &position) const
    {
      return (position(dim-1)-box_origin[dim-1]) - extents[dim-1];
    }


    template <int dim>
    Point<dim>
    Box<dim>::representative_point(const double depth) const
    {
      Assert (depth >= 0,
              ExcMessage ("Given depth must be positive or zero."));
      Assert (depth <= maximal_depth(),
              ExcMessage ("Given depth must be less than or equal to the maximal depth of this geometry."));

      // choose a point on the center axis of the domain (without topography)
      Point<dim> p = extents/2+box_origin;

      // We need a dim-1 point to get the topo value.
      Point<dim-1> surface_point;
      for (unsigned int d=0; d<dim-1; ++d)
        surface_point[d] = p[d];

      const double topo = topo_model->value(surface_point);
      p[dim-1] = extents[dim-1]+box_origin[dim-1]-depth+topo;

      return p;
    }


    template <int dim>
    double
    Box<dim>::maximal_depth() const
    {
      return extents[dim-1] + topo_model->max_topography();
    }

    template <int dim>
    std::vector<std::vector<double>>
    Box<dim>::get_surface_test() const
    {
      std::vector<std::vector<double>> surface_height;
      return surface_height;
    }

    template <int dim>
    bool
    Box<dim>::has_curved_elements() const
    {
      return false;
    }

    template <int dim>
    bool
    Box<dim>::point_is_in_domain(const Point<dim> &point) const
    {
      AssertThrow(!this->get_parameters().mesh_deformation_enabled ||
                  this->simulator_is_past_initialization() == false,
                  ExcMessage("After displacement of the free surface, this function can no longer be used to determine whether a point lies in the domain or not."));

      AssertThrow(Plugins::plugin_type_matches<const InitialTopographyModel::ZeroTopography<dim>>(this->get_initial_topography_model()),
                  ExcMessage("After adding topography, this function can no longer be used "
                             "to determine whether a point lies in the domain or not."));

      for (unsigned int d = 0; d < dim; ++d)
        if (point[d] > extents[d]+box_origin[d]+std::numeric_limits<double>::epsilon()*extents[d] ||
            point[d] < box_origin[d]-std::numeric_limits<double>::epsilon()*extents[d])
          return false;

      return true;
    }

    template <int dim>
    std::array<double,dim>
    Box<dim>::cartesian_to_natural_coordinates(const Point<dim> &position_point) const
    {
      std::array<double,dim> position_array;
      for (unsigned int i = 0; i < dim; ++i)
        position_array[i] = position_point(i);

      return position_array;
    }


    template <int dim>
    aspect::Utilities::Coordinates::CoordinateSystem
    Box<dim>::natural_coordinate_system() const
    {
      return aspect::Utilities::Coordinates::CoordinateSystem::cartesian;
    }


    template <int dim>
    Point<dim>
    Box<dim>::natural_to_cartesian_coordinates(const std::array<double,dim> &position_tensor) const
    {
      Point<dim> position_point;
      for (unsigned int i = 0; i < dim; ++i)
        position_point[i] = position_tensor[i];

      return position_point;
    }

    template <int dim>
    void
    Box<dim>::
    update_surface ()
    {
        this->get_pcout() << "   Updating surface values... " <<std::endl;
        // loop over all of the surface cells and save the elevation to a stored value.
        // This needs to be sent to 1 processor, sorted, and broadcast so that every processor knows the entire surface.
        // TODO: Is there a better time to call this in regards to mesh deformation?
        // TODO: Right now this is messy saving two variables, this can be cleaned up by
        // combining surface_x and surface_y into one variable and broadcasting them in a loop.
        std::vector<double> surface_x;
        std::vector<double> surface_y;
        std::vector<std::vector<double>> local_surface_height(2, std::vector<double>());
        const types::boundary_id relevant_boundary = this->get_geometry_model().translate_symbolic_boundary_name_to_id ("top");
        const QTrapezoid<dim-1> face_corners;
        FEFaceValues<dim> fe_face_values(this->get_mapping(),
                                          this->get_fe(),
                                          face_corners,
                                          update_quadrature_points);


        // Loop over all corners at the surface and save their X and Y positions.
        // TODO: Update this to work in 3D. Spherical?
        for (const auto &cell : this->get_dof_handler().active_cell_iterators())
          if (cell->is_locally_owned() && cell->at_boundary())
            for (const unsigned int face_no : cell->face_indices())
              if (cell->face(face_no)->at_boundary())
                {
                  if ( cell->face(face_no)->boundary_id() != relevant_boundary)
                    continue;

                  fe_face_values.reinit(cell, face_no);

                  for (unsigned int corner = 0; corner < face_corners.size(); ++corner)
                    {
                      const Point<dim> vertex = fe_face_values.quadrature_point(corner);

                      local_surface_height[0].push_back(vertex(0));   // X
                      local_surface_height[1].push_back(vertex(dim-1));  // Y, now this only works for 2D

                    }
                }

        // Combine all local_surfaces, combine and sort them, and broadcast back.
        if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          // Only push back them if the x value changes. This is done because
          // there are two corner points for each position, and we want to remove duplicates.
          // TODO: Is there an easier way to remove duplicates before getting here?
          for (unsigned int i=0; i<local_surface_height[1].size(); i++)
          {
            if(i==0)
            {
              surface_x.push_back(local_surface_height[0][i]);
              surface_y.push_back(local_surface_height[1][i]);
            }
            else
              if(local_surface_height[0][i] != local_surface_height[0][i-1])
              {
                 surface_x.push_back(local_surface_height[0][i]);
                 surface_y.push_back(local_surface_height[1][i]);
              }
          }

          for (unsigned int p=1; p<Utilities::MPI::n_mpi_processes(this->get_mpi_communicator()); ++p)
          {
            // First, find out the size of the array a process wants to send.
            MPI_Status status;
            MPI_Probe(p, 42, this->get_mpi_communicator(), &status);
            int incoming_size = 0;
            MPI_Get_count(&status, MPI_DOUBLE, &incoming_size);

            // Resize the array so it fits whatever the process sends.
            std::vector<std::vector<double>> temporary_surface(2, std::vector<double>());

            for (unsigned int i=0; i<temporary_surface.size(); ++i)
                temporary_surface[i].resize(incoming_size);

            for (unsigned int i=0; i<temporary_surface.size(); ++i)
              MPI_Recv(&temporary_surface[i][0], incoming_size, MPI_DOUBLE, p, 42, this->get_mpi_communicator(), &status);

            for (unsigned int i=0; i<temporary_surface[1].size(); ++i)
              {
                if(i==0)
                {
                  if(surface_x.size() > 0)
                  {
                    if(temporary_surface[0][i] != surface_x[surface_x.size() - 1])
                    {
                      surface_x.push_back(temporary_surface[0][i]);
                      surface_y.push_back(temporary_surface[1][i]);
                    }
                  else
                    {
                      surface_x.push_back(temporary_surface[0][i]);
                      surface_y.push_back(temporary_surface[1][i]);
                    }
                  }
                }
                else
                {
                  if(temporary_surface[0][i] != temporary_surface[0][i-1])
                  {
                    surface_x.push_back(temporary_surface[0][i]);
                    surface_y.push_back(temporary_surface[1][i]);
                  }
                }
              }
          }

            double vector_size = surface_x.size();
            MPI_Bcast(&vector_size, 1, MPI_DOUBLE, 0, this->get_mpi_communicator());
            MPI_Bcast(&surface_x[0], vector_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
            MPI_Bcast(&surface_y[0], vector_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }
        else
        {
          for (unsigned int i=0; i<local_surface_height.size(); i++)
              MPI_Ssend(&local_surface_height[i][0], local_surface_height[1].size(), MPI_DOUBLE, 0, 42, this->get_mpi_communicator());

          double vector_size = 0;
         //std::vector<std::vector<double>> surface_height(2, std::vector<double>());
          MPI_Bcast(&vector_size, 1, MPI_DOUBLE, 0, this->get_mpi_communicator());
          surface_x.resize(vector_size);
          surface_y.resize(vector_size);

          MPI_Bcast(&surface_x[0], vector_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
          MPI_Bcast(&surface_y[0], vector_size, MPI_DOUBLE, 0, this->get_mpi_communicator());
        }

        //surface_height[0] = surface_x;
        //surface_height[1] = surface_y;

        surface_xx = surface_x;
        surface_yy = surface_y;
    }

    template <int dim>
    double
    Box<dim>::depth_including_mesh_deformation(const Point<dim> &position) const
    {
          // Compute the depth including free surface position, this only works in 2D
          // The surface points do not always match position points, so we linearly
          // Interpolate between the two nearest surface points.
          // This function is sometimes called before the surface has been saved,
          // so during timestep zero we call the normal depth function.
          double depth_from_surface;
          if(this->get_timestep_number() > 0)
          {
            double height = 0;
            double x1 = 0; double x2 = 0; double y1 = 0; double y2 = 0;

            // Loop through all surface points and find the nearest ones to the given position.
            for(unsigned int s=0; s<surface_xx.size();++s)
            {
              x2 = surface_xx[s];
              y2 = surface_yy[s];

              if(x2 > position[0])
                break;

              // If we aren't greater, then increase point1.
              x1 = x2; y1 = y2;
            }

            // If we are at the edge of the model domain, use that y value. Otherwise interpolate.
            if(position[0] >= surface_xx[surface_xx.size() - 1])
              height = surface_yy[surface_yy.size() - 1];
            else if (position[0] <= surface_xx[0])
              height = surface_yy[0];
            else
              height = ((y2 - y1)/(x2-x1))*(position[0] - x1) + y1;

            depth_from_surface = height - position[dim-1];
          }
          else
            depth_from_surface = depth(position);

          return std::min (std::max (depth_from_surface, 0.), maximal_depth());
    }


    template <int dim>
    void
    Box<dim>::
    declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          prm.declare_entry ("X extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in x-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Y extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in y-direction. Units: \\si{\\meter}.");
          prm.declare_entry ("Z extent", "1.",
                             Patterns::Double (0.),
                             "Extent of the box in z-direction. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("Box origin X coordinate", "0.",
                             Patterns::Double (),
                             "X coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Y coordinate", "0.",
                             Patterns::Double (),
                             "Y coordinate of box origin. Units: \\si{\\meter}.");
          prm.declare_entry ("Box origin Z coordinate", "0.",
                             Patterns::Double (),
                             "Z coordinate of box origin. This value is ignored "
                             "if the simulation is in 2d. Units: \\si{\\meter}.");

          prm.declare_entry ("X repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in X direction.");
          prm.declare_entry ("Y repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Y direction.");
          prm.declare_entry ("Z repetitions", "1",
                             Patterns::Integer (1),
                             "Number of cells in Z direction.");

          prm.declare_entry ("X periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in X direction");
          prm.declare_entry ("Y periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Y direction");
          prm.declare_entry ("Z periodic", "false",
                             Patterns::Bool (),
                             "Whether the box should be periodic in Z direction");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    Box<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Box");
        {
          box_origin[0] = prm.get_double ("Box origin X coordinate");
          extents[0] = prm.get_double ("X extent");
          periodic[0] = prm.get_bool ("X periodic");
          repetitions[0] = prm.get_integer ("X repetitions");

          if (dim >= 2)
            {
              box_origin[1] = prm.get_double ("Box origin Y coordinate");
              extents[1] = prm.get_double ("Y extent");
              periodic[1] = prm.get_bool ("Y periodic");
              repetitions[1] = prm.get_integer ("Y repetitions");
            }

          if (dim >= 3)
            {
              // Use dim-1 instead of 2 to avoid compiler warning in 2d:
              box_origin[dim-1] = prm.get_double ("Box origin Z coordinate");
              extents[dim-1] = prm.get_double ("Z extent");
              periodic[dim-1] = prm.get_bool ("Z periodic");
              repetitions[dim-1] = prm.get_integer ("Z repetitions");
            }
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
  namespace GeometryModel
  {
    ASPECT_REGISTER_GEOMETRY_MODEL(Box,
                                   "box",
                                   "A box geometry parallel to the coordinate directions. "
                                   "The extent of the box in each coordinate direction "
                                   "is set in the parameter file. The box geometry labels its "
                                   "2*dim sides as follows: in 2d, boundary indicators 0 through 3 "
                                   "denote the left, right, bottom and top boundaries; in 3d, boundary "
                                   "indicators 0 through 5 indicate left, right, front, back, bottom "
                                   "and top boundaries (see also the documentation of the deal.II class "
                                   "``ReferenceCell''). You can also use symbolic names ``left'', ``right'', "
                                   "etc., to refer to these boundaries in input files. "
                                   "It is also possible to add initial topography to the box model. Note however that "
                                   "this is done after the last initial adaptive refinement cycle. "
                                   "Also, initial topography is supposed to be small, as it is not taken into account "
                                   "when depth or a representative point is computed. ")


  }
}
