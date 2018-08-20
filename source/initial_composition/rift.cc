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
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
 */


#include <aspect/initial_composition/rift.h>
#include <aspect/postprocess/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/ellipsoidal_chunk.h>
#include <aspect/material_model/visco_plastic.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    Rift<dim>::Rift ()
    {}

    template <int dim>
    void
    Rift<dim>::initialize ()
    {
      AssertThrow(dynamic_cast<const MaterialModel::ViscoPlastic<dim> *>(&this->get_material_model()) != NULL,
                  ExcMessage("This initial condition only makes sense in combination with the visco_plastic material model."));

      // From shear_bands.cc
      Point<dim> extents_min, extents_max;
      TableIndices<dim> size_idx;
      for (unsigned int d=0; d<dim; ++d)
        size_idx[d] = grid_intervals[d]+1;

      Table<dim,double> white_noise;
      white_noise.TableBase<dim,double>::reinit(size_idx);
      std_cxx1x::array<std::pair<double,double>,dim> grid_extents;

      if (cartesian_domain)
        {
          const GeometryModel::Box<dim> *geometry_model
            = dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model());

          // Min and max of each direction (m)
          extents_max = geometry_model->get_extents();
        }
      else if (const GeometryModel::Chunk<dim> *geometry_model
               = dynamic_cast<const GeometryModel::Chunk<dim> *>(&this->get_geometry_model()))
        {
          // Min and max radius
          extents_max[0] = geometry_model->outer_radius();
          extents_min[0] = std::max(extents_max[0]-strain_depth-5.*strain_halfwidth,geometry_model->inner_radius());

          // Min and max longitude (already in radians)
          extents_min[1] = geometry_model->east_longitude();
          extents_max[1] = geometry_model->west_longitude();

          // Min and max latitude (already in radians), convert to colatitude
          if (dim == 3)
            {
              extents_max[dim-1] = 0.5 * numbers::PI - geometry_model->south_latitude();
              extents_min[dim-1] = 0.5 * numbers::PI - geometry_model->north_latitude();
            }
        }
      else if (const GeometryModel::EllipsoidalChunk<dim> *geometry_model
               = dynamic_cast<const GeometryModel::EllipsoidalChunk<dim> *>(&this->get_geometry_model()))
        {
          // Check that the model is not elliptical
          AssertThrow(geometry_model->get_eccentricity() == 0.0, ExcMessage("This boundary velocity plugin cannot be used with a non-zero eccentricity. "));

          // Min and max radius
          extents_max[0] = geometry_model->get_semi_major_axis_a();
          extents_min[0] = std::max(extents_max[0]-geometry_model->maximal_depth(), extents_max[0]-strain_depth-5.*strain_halfwidth);
          // Assume chunk outlines are lat/lon parallel
          std::vector<Point<2> > corners = geometry_model->get_corners();
          // Convert to radians, lon, colat
          extents_min[1] = corners[1][0] * numbers::PI / 180.;
          extents_max[1] = corners[0][0] * numbers::PI / 180.;
          extents_min[dim-1] = 0.5 * numbers::PI - corners[0][1] * numbers::PI / 180.;
          extents_max[dim-1] = 0.5 * numbers::PI - corners[2][1] * numbers::PI / 180.;
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("This initial condition only works with the box or (ellipsoidal) chunk geometry model."));
        }

      for (unsigned int d=0; d<dim; ++d)
        {
          grid_extents[d].first=extents_min[d];
          grid_extents[d].second=extents_max[d];
        }

      // use a fixed number as seed for random generator
      // this is important if we run the code on more than 1 processor
      std::srand(seed);

      TableIndices<dim> idx;

      for (unsigned int i=0; i<white_noise.size()[0]; ++i)
        {
          idx[0] = i;
          for (unsigned int j=0; j<white_noise.size()[1]; ++j)
            {
              idx[1] = j;
              if (dim == 3)
                {
                  for (unsigned int k=0; k<white_noise.size()[dim-1]; ++k)
                    {
                      idx[dim-1] = k;
                      // std::rand will give a value between zero and RAND_MAX (usually INT_MAX).
                      // The modulus of this value and 10000, gives a value between 0 and 10000-1.
                      // Subsequently dividing by 5000.0 will give value between 0 and 2 (excluding 2).
                      // Subtracting 1 will give a range [-1,1)
                      // Because we want values [0,1), we change our white noise computation to:
                      white_noise(idx) = ((std::rand() % 10000) / 10000.0);
                    }
                }
              else
                white_noise(idx) = ((std::rand() % 10000) / 10000.0);
            }
        }

      interpolate_noise = new Functions::InterpolatedUniformGridData<dim> (grid_extents,
                                                                           grid_intervals,
                                                                           white_noise);
    }

    template <int dim>
    double
    Rift<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      // If n_comp does not represent the strain field,
      // return 0 right away.
      if (n_comp != strain_composition_number)
        return 0.0;

      // Initiate distance with large value
      double distance_to_rift_axis = 1e23;
      double temp_distance = 0;

      // For spherical geometries we need to reorder the coordinates
      Point<dim> natural_coords = position;

      // Loop over all line segments
      for (unsigned int i_segments = 0; i_segments < point_list.size(); ++i_segments)
        {
          if (cartesian_domain)
            {
              if (dim == 2)
                temp_distance = std::abs(natural_coords[0]-point_list[i_segments][0][0]);
              else
                {
                  // Get the surface coordinates by dropping the last coordinate
                  const Point<2> surface_position = Point<2>(natural_coords[0],natural_coords[1]);
                  temp_distance = std::abs(Utilities::distance_to_line(point_list[i_segments], surface_position));
                }
            }
          // chunk (spherical) geometries
          else
            {
              // spherical coordinates in radius [m], lon [rad], colat [rad] format
              const std_cxx11::array<double,dim> spherical_point = Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
              natural_coords[0] = spherical_point[0];
              Point<2> surface_position;
              for (unsigned int d=0; d<dim-1; ++d)
                {
                  surface_position[d] = spherical_point[d+1];
                  natural_coords[d+1] = spherical_point[d+1];
                }

              temp_distance = (dim == 2) ? std::abs(surface_position[0]-point_list[i_segments][0][0]) : Utilities::distance_to_line(point_list[i_segments], surface_position);
            }

          // Get the minimum distance
          distance_to_rift_axis = std::min(distance_to_rift_axis, temp_distance);
        }

      // Smoothing of noise with depth
      const double depth_smoothing = 0.5 * (1.0 - std::tanh((this->get_geometry_model().depth(position) - strain_depth) / strain_halfwidth));
      // Smoothing of noise with lateral distance to the rift axis
      const double noise_amplitude = A * std::exp((-std::pow(distance_to_rift_axis,2)/(2.0*std::pow(sigma,2)))) * depth_smoothing;
      // Add randomness
      return noise_amplitude * interpolate_noise->value(natural_coords);
    }

    template <int dim>
    void
    Rift<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Rift");
        {
          prm.declare_entry ("Random number generator seed", "0",
                             Patterns::Double (0),
                             "The value of the seed used for the random number generator. "
                             "Units: none.");
          prm.declare_entry ("Standard deviation of Gaussian noise amplitude distribution", "20000",
                             Patterns::Double (0),
                             "The standard deviation of the Gaussian distribution of the amplitude of the strain noise. "
                             "Note that this parameter is taken to be the same for all rift segments. "
                             "Units: $m$ or degrees.");
          prm.declare_entry ("Maximum amplitude of Gaussian noise amplitude distribution", "0.2",
                             Patterns::Double (0),
                             "The amplitude of the Gaussian distribution of the amplitude of the strain noise. "
                             "Note that this parameter is taken to be the same for all rift segments. "
                             "Units: none.");
          prm.declare_entry ("Depth around which Gaussian noise is smoothed out", "40000",
                             Patterns::Double (0),
                             "The depth around which smoothing out of the strain noise starts with a hyperbolic tangent. "
                             "Note that this parameter is taken to be the same for all rift segments. "
                             "Units: $m$.");
          prm.declare_entry ("Halfwidth with which Gaussian noise is smoothed out in depth", "40000",
                             Patterns::Double (0),
                             "The halfwidth with which smoothing out of the strain noise is done with a hyperbolic tangent. "
                             "Note that this parameter is taken to be the same for all rift segments. "
                             "Units: $m$.");
          prm.declare_entry ("Grid intervals for noise X or radius", "25",
                             Patterns::Integer (0),
                             "Grid intervals in X (cartesian domain) or radial (spherical) direction for the white noise "
                             "added to the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
          prm.declare_entry ("Grid intervals for noise Y or longitude", "25",
                             Patterns::Integer (0),
                             "Grid intervals in Y (cartesian domain) or longitude (spherical) direction for the white noise "
                             "added to the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
          prm.declare_entry ("Grid intervals for noise Z or latitude", "25",
                             Patterns::Integer (0),
                             "Grid intervals in Z (cartesian domain) or latitude (spherical) direction for the white noise "
                             "added to the initial background porosity that will then be interpolated "
                             "to the model grid. "
                             "Units: none.");
          prm.declare_entry ("Rift axis line segments",
                             "",
                             Patterns::Anything(),
                             "Set the line segments that represent the rift axis. In 3d each segment is made up of "
                             "two points that represent horizontal coordinates (x,y) or (lon,lat). "
                             "The exact format for the point list describing the segments is "
                             "\"x1,y1>x2,y2;x2,y2>x3,y3;x4,y4>x5,y5\". In 2d, a segment is made up by 1 horizontal "
                             "x or longitude coordinate: \"x1;x2;x3\". Note that the segments can be connected "
                             "or isolated. The units of the coordinates are "
                             "dependent on the geometry model. In the box model they are in meters, in the "
                             "chunks they are in degrees.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Rift<dim>::parse_parameters (ParameterHandler &prm)
    {
      // Check that there is a compositional field called strain and retrieve its index
      
      if(this->introspection().compositional_name_exists("plastic_strain"))
                  strain_composition_number = this->introspection().compositional_index_for_name("plastic_strain");
      else if(this->introspection().compositional_name_exists("viscous_strain"))
                  strain_composition_number = this->introspection().compositional_index_for_name("viscous_strain");
      else if(this->introspection().compositional_name_exists("total_strain"))
                  strain_composition_number = this->introspection().compositional_index_for_name("total_strain");
      else
          AssertThrow(false, ExcMessage("This plugin requires a compositional strain field (plastic_strain, viscous_strain, or total_strain). "));

      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Rift");
        sigma                = prm.get_double ("Standard deviation of Gaussian noise amplitude distribution");
        A                    = prm.get_double ("Maximum amplitude of Gaussian noise amplitude distribution");
        seed                 = prm.get_double ("Random number generator seed");
        strain_depth         = prm.get_double ("Depth around which Gaussian noise is smoothed out");
        strain_halfwidth     = prm.get_double ("Halfwidth with which Gaussian noise is smoothed out in depth");
        grid_intervals[0]    = prm.get_integer ("Grid intervals for noise X or radius");
        grid_intervals[1]    = prm.get_integer ("Grid intervals for noise Y or longitude");
        if (dim == 3)
          grid_intervals[2]    = prm.get_integer ("Grid intervals for noise Z or latitude");

        // Read in the string of rift segments
        const std::string temp_all_segments = prm.get("Rift axis line segments");
        // Split the string into segment strings
        const std::vector<std::string> temp_segments = Utilities::split_string_list(temp_all_segments,';');
        // The number of segments, each consisting of a begin and an end point in 3d and one point in 3d
        const unsigned int n_temp_segments = temp_segments.size();
        point_list.resize(n_temp_segments);

        // Default is true, but in case we use a (ellipsoidal) chunk domain, set to false
        if (dynamic_cast<const GeometryModel::Box<dim> *>(&this->get_geometry_model()) == NULL)
          cartesian_domain = false;

        // Loop over the segments to extract the points
        for (unsigned int i_segment = 0; i_segment < n_temp_segments; i_segment++)
          {
            // In 3d a line segment consists of 2 points,
            // in 2d of only 1 (ridge axis orthogonal to x and y).
            // Also, a 3d point has 2 coordinates (x and y),
            // a 2d point only 1 (x).
            //point_list[i_segment].resize(dim-1);

            if (dim == 3)
              {
                const std::vector<std::string> temp_segment = Utilities::split_string_list(temp_segments[i_segment],'>');
                Assert(temp_segment.size() == 2,ExcMessage ("The given coordinate '" + temp_segment[i_segment] + "' is not correct. "
                                                            "It should only contain 2 parts: "
                                                            "the two points of the segment, separated by a '>'."));

                // Loop over the 2 points of each segment
                for (unsigned int i_points = 0; i_points < dim-1; i_points++)
                  {
                    std::vector<double> temp_point = Utilities::string_to_double(Utilities::split_string_list(temp_segment[i_points],','));
                    Assert(temp_point.size() == 2,ExcMessage ("The given coordinates of segment '" + temp_segment[i_points] + "' are not correct. "
                                                              "It should only contain 2 parts: "
                                                              "the x and y coordinates of the segment begin/end point, separated by a ','."));

                    if (!cartesian_domain)
                      {
                        // convert degrees to radians for (ellipsoidal) chunks
                        // longitude
                        temp_point[0] *= numbers::PI/180.;
                        // and convert latitude to colatitude
                        temp_point[1] = 0.5 * numbers::PI - temp_point[1] * numbers::PI / 180.;
                      }

                    // Add the point to the list of points for this segment
                    point_list[i_segment][i_points] = (Point<2>(temp_point[0], temp_point[1]));
                  }
              }
            else
              {
                // Add the point to the list of points for this segment
                // As we're in 2d all segments correspond to 1 point consisting of 1 coordinate
                const double temp_point = Utilities::string_to_double(temp_segments[i_segment]);
                point_list[i_segment][0] = (Point<2>(temp_point, temp_point));
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
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Rift,
                                              "rift",
                                              "Specify the strain initial compositional field value based on the distance to a list of line segments "
                                              "and the user-defined Gaussian distribution around these segments, combined with random noise.")
  }
}
