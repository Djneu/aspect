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

#ifndef _aspect_mesh_deformation_fastscape_h
#define _aspect_mesh_deformation_fastscape_h

#include <aspect/global.h>

#ifdef ASPECT_WITH_FASTSCAPE

#include <aspect/mesh_deformation/interface.h>

namespace aspect
{
  using namespace dealii;

  namespace MeshDeformation
  {
    /**
     * Define FastScape functions as C functions. Must use exact same function/variable name
     * and type as used in FastScape. All function names must be made lowercase, and an
     * underscore added at the end. Types must be defined as pointers, and sent to
     * FastScape as a reference. Additional functions are available within FastScape,
     * see https://fastscape.org/fastscapelib-fortran/ for a list of all functions.
     */
#ifdef __cplusplus
    extern"C" {
#endif
    /**
     * Function to initialize FastScape.
     */
    void fastscape_init_();

    /**
     * Set the x and y extent of the FastScape model.
     */
    void fastscape_set_xl_yl_(const double *xxl,const double *yyl);

    /**
     * Set number of grid points in x (nx) and y (ny)
     */
    void fastscape_set_nx_ny_(const int *nnx, const int *nny);

    /**
     * Allocate memory, must be called after set nx/ny.
     */
    void fastscape_setup_();

    /**
     * Set FastScape boundary conditions.
     */
    void fastscape_set_bc_(const int *jbc);

    /**
     * Set FastScape timestep. This will vary based on the ASPECT timestep.
     */
    void fastscape_set_dt_(double *dtt);

    /**
     * Initialize FastScape topography.
     */
    void fastscape_init_h_(double *hp);

    /**
     * Set FastScape erosional parameters on land. These parameters will apply to the stream power law (SPL)
     * and hillslope diffusion for basement and sediment. This can be set between timesteps.
     */
    void fastscape_set_erosional_parameters_(double *kkf,const double *kkfsed,const double *mm,const double *nnn,
                                             double *kkd,const double *kkdsed,const double *gg1,const double *gg2,const double *pp);

    /**
     * Set FastScape marine erosional parameters. This can be set between timesteps.
     */
    void fastscape_set_marine_parameters_(const double *sl, const double *p1, const double *p2, const double *z1,
                                          const double *z2, const double *r, const double *l, const double *kds1, const double *kds2);

    /**
     * Set advection velocities for FastScape. This can be set between timesteps.
     */
    void fastscape_set_v_(double *ux, double *uy);

    /**
     * Set FastScape uplift rate. This can be set between timesteps.
     */
    void fastscape_set_u_(double *up);

    /**
     * Set FastScape topography. This can be set between timesteps.
     */
    void fastscape_set_h_(double *hp);

    /**
    * Set FastScape basement. This can be set between timesteps. Sediment within FastScape
    * is considered as the difference between the topography and basement, though this may differ
    * from sediment as seen in ASPECT because the FastScape basement only takes the surface
    * velocities into consideration.
    */
    void fastscape_set_basement_(double *b);

    /**
     * Get the current FastScape step iteration. This is set to zero when restarting FastScape.
     */
    void fastscape_get_step_(int *sstep);

    /**
     * Run FastScape for a single FastScape timestep.
     */
    void fastscape_execute_step_();

    /**
     * Create a .VTK file for the FastScape surface within the VTK folder of the
     * ASPECT output folder.
     */
    void fastscape_named_vtk_(double *fp, const double *vexp, int *astep, const char *c, int *length);

    /**
     * Copy the current FastScape topography.
     */
    void fastscape_copy_h_(double *hp);

    /**
     * Copy the current FastScape basement.
     */
    void fastscape_copy_basement_(double *b);

    /**
     * Initialize the straigraphy component of FastScape.
     */
    void fastscape_strati_(const int *nstepp, const int *nreflectorp, int *nfreqp, const double *vexp);

    /**
     * Copy the current FastScape slopes.
     */
    void fastscape_copy_slope_(double *slopep);

    /**
     * Destroy FastScape.
     */
    void fastscape_destroy_();
#ifdef  __cplusplus
  }
#endif

  /**
   * A plugin that utilizes the landscape evolution code FastScape
   * to deform the ASPECT boundary through advection, uplift,
   * hillslope diffusion, the stream power law, sediment deposition
   * and marine diffusion.
   */
  template<int dim>
  class FastScape : public Interface<dim>, public SimulatorAccess<dim>
  {
    public:
      /**
       * Initialize variables for FastScape.
       */
      virtual void initialize ();

      /**
        * A function that creates constraints for the velocity of certain mesh
        * vertices (e.g. the surface vertices) for a specific boundary.
        * The calling class will respect
        * these constraints when computing the new vertex positions.
        */
      virtual
      void
      compute_velocity_constraints_on_boundary(const DoFHandler<dim> &mesh_deformation_dof_handler,
                                               AffineConstraints<double> &mesh_velocity_constraints,
                                               const std::set<types::boundary_id> &boundary_id) const;

      /**
       * Declare parameters for the FastScape plugin.
       */
      static
      void declare_parameters (ParameterHandler &prm);

      /**
       * Parse parameters for the FastScape plugin.
       */
      void parse_parameters (ParameterHandler &prm);

      /**
       * Function used to set the FastScape ghost nodes. FastScape boundaries are
       * not uplifted or periodic for advection and diffusion. By using a layer
       * of extra nodes in the FastScape model, we can avoid seeing FastScape boundary
       * effects within ASPECT. Similalrly, we can use these nodes to have fully
       * periodic boundaries, where we check the flow direction and update the FastScape
       * ghost nodes, and the nodes one layer inward (boundary nodes in ASPECT) to match
       * the parameters on the other side (vx, vy, vz, h). This is done every ASPECT timestep
       * before running FastScape.
       * TODO: Should this be in private?
       */
      void set_ghost_nodes(double *h, double *vx, double *vy, double *vz) const;

      /**
       * Function to fill the Fastscape arrays (height and velocities) with the data received from ASPECT in the correct index order.
       */
      void fill_fastscape_arrays(double *h, double *kd, double *kf, double *vx, double *vy, double *vz, std::vector<std::vector<double>> temporary_variables) const;

      /**
       * Function to get the ASPECT topography and velocities at the surface, and an index for transferring these to FastScape.
       */
      std::vector<std::vector<double>> get_aspect_values() const;

      /**
       * Function to intialize or restart FastScape
       */
      void initialize_fastscape(double *h, double *b, double *kd, double *kf) const;

      /**
       * Execute FastScape
       */
      void execute_fastscape(double *h, double *kd, int istep) const;

      /**
       * Function to apply orographic controls to the FastScape model.
       */
      void apply_orographic_controls(double *h, double *kf, double *kd) const;

      /**
       * Fill velocity data table to be interpolated back onto the ASPECT mesh.
       */
      Table<dim,double> fill_data_table(std::vector<double> values, TableIndices<dim> size_idx, int nx, int ny) const;

    private:
      /**
       * First attempt at the number of steps to run FastScape for every ASPECT timestep,
       * where the FastScape timestep is determined by (ASPECT_timestep/nstep).
       */
      int nstep;

      /**
       * Maximum timestep allowed for FastScape, if (ASPECT_timestep/nsteps) exceeds this, nsteps is doubled.
       */
      double maximum_fastscape_timestep;

      /**
       * Check whether FastScape needs to be restarted. This is used as
       * a mutable bool because we determine whether the model is being resumed in
       * initialize(), and then after reinitializing FastScape we change it to zero
       * so it does not initialize FastScape again in future timesteps.
       * TODO: There is probably a better way to do this, and restarts should be rolled into
       * the general ASPECT restart.
       */
      mutable bool restart;

      /**
       * FastScape resets to step 0 on restart, this keeps the step number from the previous run.
       * As we only want to read this value in from the restart file once, we initialize it in
       * initialize() and only change it during the first timestep where FastScape is initialized.
       */
      mutable int restart_step;

      /**
       * ASPECT end time to check if we should destroy FastScape.
       */
      double end_time;

      /**
       * FastScape cell size in X, dx should always equal dy.
       */
      double dx;

      /**
       * FastScape cell size in Y, dy should always equal dx.
       */
      double dy;

      /**
       * FastScape X extent (ASPECT X extent + 2*dx for ghost nodes).
       */
      double x_extent;

      /**
       * Fastscape Y extent (ASPECT Y extent + 2*dy for ghost nodes).
       */
      double y_extent;

      /**
       * User set FastScape Y extent for a 2D ASPECT model.
       */
      double y_extent_2d;

      /**
       * Number of x points in FastScape array.
       */
      int nx;

      /**
       * Number of y points in FastScape array.
       */
      int ny;

      /**
       * Size of the FastScape array (nx*ny).
       */
      int array_size;

      /**
       * Vertical exaggeration in FastScape visualization.
       */
      double vexp;

      /**
       * How many levels FastScape should be refined above the maximum ASPECT surface resolution.
       */
      unsigned int additional_refinement;

      /**
       * Maximum expected refinement level at ASPECT's surface.
       * This and resolution_difference are required to properly transfer node data from
       * ASPECT to FastScape.
       */
      int maximum_surface_refinement_level;

      /**
       * Difference in refinement levels expected at the ASPECT surface,
       * where this would be set to 2 if 3 refinement leves are set at the surface.
       * This and surface_resolution are required to properly transfer node data from
       * ASPECT to FastScape.
       *
       * TODO: Should this be kept this way, or make it so the input is the expected levels
       * of refinement at the surface, and we can subtract one within the code? Also,
       * it would be good to find a way to check these are correct, because they are a
       * common source of errors.
       */
      int surface_refinement_difference;

      /**
       * If set to false, the FastScape surface is averaged along Y and returned
       * to ASPECT. If set to true, the center slice of the FastScape model is
       * returned to ASPECT.
       *
       * TODO: Do we average the ghost nodes or remove them? Need to double check.
       */
      bool center_slice;

      /**
       * Seed number for initial topography noise in FastScape.
       */
      int fs_seed;

      /**
       * Variable to hold ASPECT domain extents.
       */
      std::array<std::pair<double,double>,dim> grid_extent;

      /**
       * Table for interpolating FastScape surface velocities back to ASPECT.
       */
      std::array< unsigned int, dim > table_intervals;

      /**
       * Check whether FastScape is initialized at the surface.
       */
      //std::map<types::boundary_id, std::vector<std::string>> mesh_deformation_boundary_indicators_map;

      /**
       * Whether or not to use the ghost nodes.
       */
      bool use_ghost_nodes;

      /**
       * Magnitude (m) of the initial noise applied to FastScape.
       * Applied as either a + or - value to the topography
       * such that the total difference can be up to 2*noise_h.
       */
      int noise_h;

      /**
       * Sediment rain in m/yr, added as a flat increase to the FastScape surface
       * every ASPECT timestep before running FastScape.
       */
      std::vector<double> sediment_rain_rates;

      /**
       * Time at which each interval of sediment_rain_rates is active. Should contain
       * one less value than sediment_rain_rates, assuming that sediment_rain_rates[0]
       * is applied from model time 0 until sediment_rain_times[0].
       */
      std::vector<double> sediment_rain_times;

      // FastScape boundary condition variables //
      /**
       * FastScape bottom boundary condition that determines topography at the FastScape bottom boundary.
       * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
       * reflective boundary. When two opposing boundaries are relfective (e.g., top and bottom are both zero), then the boundaries
       * become cyclic.
       */
      unsigned int bottom;

      /**
       * FastScape top boundary condition that determines topography at the FastScape top boundary.
       * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
       * reflective boundary. When two opposing boundaries are relfective (e.g., top and bottom are both zero), then the boundaries
       * become cyclic.
       */
      unsigned int top;

      /**
       * FastScape right boundary condition that determines topography at the FastScape right boundary.
       * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
       * reflective boundary. When two opposing boundaries are relfective (e.g., left and right are both zero), then the boundaries
       * become cyclic.
       */
      unsigned int right;

      /**
       * FastScape left boundary condition that determines topography at the FastScape left boundary.
       * Where 1 represents a fixed height boundary (though this can still be uplifted through uplift velocities), and 0 a
       * reflective boundary. When two opposing boundaries are relfective (e.g., left and right are both zero), then the boundaries
       * become cyclic.
       */
      unsigned int left;

      /**
       * Paramters that set the fastscape boundaries periodic even though the ghost nodes are set 'fixed'
       */
      bool topbottom_ghost_nodes_periodic;
      bool leftright_ghost_nodes_periodic;

      /**
       * Integer that holds the full boundary conditions sent to FastScape (e.g., 1111).
       */
      int bc;

      /**
       * Prescribed flux per unit length into the model through the bottom boundary (m^2/yr).
       */
      double bottom_flux;

      /**
       * Prescribed flux per unit length into the model through the top boundary (m^2/yr).
       */
      double top_flux;

      /**
       * Prescribed flux per unit length into the model through the right boundary (m^2/yr).
       */
      double right_flux;

      /**
       * Prescribed flux per unit length into the model through the left boundary (m^2/yr).
       */
      double left_flux;

      // FastScape erosional parameters //
      /**
       * Drainage area exponent for the stream power law (m variable in FastScape surface equation).
       */
      double m;

      /**
       * Slope exponent for the steam power law (n variable in FastScape surface equation).
       */
      double n;

      /**
       * Slope exponent for multi-direction flow, where 0 is uniform, and 10 is steepest descent. (-1 varies with slope)
       * (p variable in FastScape surface equation).
       */
      double p;

      /**
       * Bedrock deposition coefficient. Higher values deposit more sediment
       * inside the domain. (G variable in FastScape surface equation).
       */
      double g;

      /**
       * Sediment deposition coefficient. Higher values deposit more sediment inside the domain.
       * When set to -1 this is identical to the bedrock value.
       * (G variable in FastScape surface equation applied to sediment).
       */
      double gsed;

      /**
       * Bedrock river incision rate for the stream power law
       * (meters^(1-2m)/yr, kf variable in FastScape surface equation).
       */
      double kff;

      /**
       * Sediment river incision rate for the stream power law (meters^(1-2m)/yr).
       * When set to -1 this is identical to the bedrock value.
       * (kf variable in FastScape surface equation applied to sediment).
       */
      double kfsed;

      /**
       * Bedrock transport coefficient for hillslope diffusion (m^2/yr, kd in FastScape surface equation).
       */
      double kdd;

      /**
       * Bedrock transport coefficient for hillslope diffusion (m^2/yr). When set to -1 this is
       * identical to the bedrock value.
       * (kd in FastScape surface equation applied to sediment).
       */
      double kdsed;

      // FastScape marine parameters //
      /**
       * Fastscape sea level (m), set relative to the ASPECT surface where
       * a sea level of zero will represent the maximum Y (2D) or Z (3D) extent
       * inside ASPECT.
       */
      double sl;

      // Parameters to set an extra erosional base level
      // on the ghost nodes that differs from sea level.
      bool use_extra_base_level;
      double h_extra_base_level;

      /**
       * Surface porosity for sand.
       */
      double p1;

      /**
       * Surface porosity for shale.
       */
      double p2;

      /**
       * Sands e-folding depth for exponential porosity law (m).
       */
      double z1;

      /**
       * Shales e-folding depth for exponential porosity law (m).
       */
      double z2;

      /**
       * Sand-shale ratio
       */
      double r;

      /**
       * Averaging depth/thickness for sand-shale equation (m).
       */
      double l;

      /**
       * Sand marine transport coefficient (marine diffusion, m^2/yr).
       */
      double kds1;

      /**
       * Shale marine transport coefficient (marine diffusion, m^2/yr).
       */
      double kds2;

      /**
       * Flag to use the marine component of FastScape.
       */
      bool use_marine;

      // Orographic parameters //
      /**
       * Set a flat height (m) after which the flat_erosional_factor
       * is applied to the kf and kd terms.
       */
      int flat_elevation;

      /**
       * Set the height (m) after which the model will track the ridge line
       * and based on the wind direction will apply the wind_barrier_erosional_factor
       * to the kf and kd terms.
       */
      int wind_barrier_elevation;

      /**
       * Wind direction for wind_barrier_erosional_factor.
       */
      int wind_direction;

      /**
       * Factor to multiply kf and kd by depending on them
       * flat_elevation.
       */
      double flat_erosional_factor;

      /**
       * Factor to multiply kf and kd by depending on them
       * wind_barrier_elevation and wind direction.
       */
      double wind_barrier_erosional_factor;

      /**
       * Flag to stack both orographic controls.
       */
      bool stack_controls;

      /**
       * Flag to use orographic controls.
       */
      bool use_orographic_controls;

      /**
       * Flag to use stratigraphic component of FastScape.
       */
      bool use_stratigraphy;

      /**
       *
       */
      int nstepp;

      /**
       *
       */
      int nreflectorp;

      /**
       * Flag for having FastScape advect/uplift the surface. If the free surface is used
       * in conjuction with FastScape, this can be set to false, and  then FastScape will only
       * apply erosion/deposition to the surface and not advect or uplift it.
       */
      bool use_velocities;

      /**
       * Precision value for how close a ASPECT node must be to the FastScape node
       * for the value to be transferred. This is only necessary if use_v is set to 0
       * and the free surface is used to advect the surface with a normal projection.
       */
      double precision;

      /**
        * Interval between the generation of graphical output. This parameter
        * is read from the input file and consequently is not part of the
        * state that needs to be saved and restored.
        */
      double output_interval;

      /**
        * A time (in seconds) at which the last graphical output was supposed
        * to be produced. Used to check for the next necessary output time.
        */
      mutable double last_output_time;
  };
}
}

#endif
#endif
