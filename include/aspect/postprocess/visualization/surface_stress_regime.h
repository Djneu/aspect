/*
  Copyright (C) 2016 - 2017 by the authors of the ASPECT code.
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


#ifndef _aspect_postprocess_visualization_surface_stress_regime_h
#define _aspect_postprocess_visualization_surface_stress_regime_h

#include <aspect/postprocess/visualization.h>
#include <aspect/simulator_access.h>

#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/base/symmetric_tensor.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      enum struct SurfaceSymmetricTensorEigenvectorMethod
      {
        /**
         * A hybrid approach that preferentially uses the characteristic equation to
         * compute eigenvalues and an analytical approach based on the cross-product
         * for the eigenvectors. If the computations are deemed too inaccurate then
         * the method falls back to ql_implicit_shifts.
         *
         * This method potentially offers the quickest computation if the pathological
         * case is not encountered.
         */
        hybrid,
        /**
         * The iterative QL algorithm with implicit shifts applied after
         * tridiagonalization of the tensor using the householder method.
         *
         * This method offers a compromise between speed of computation and its
         * robustness. This method is particularly useful when the elements
         * of $T$ have greatly varying magnitudes, which would typically lead to a
         * loss of accuracy when computing the smaller eigenvalues.
         */
        ql_implicit_shifts,
        /**
         * The iterative Jacobi algorithm.
         *
         * This method offers is the most robust of the available options, with
         * reliable results obtained for even the most pathological cases. It is,
         * however, the slowest algorithm of all of those implemented.
         */
        jacobi
      };

      namespace Surfaceinternal
      {
        // I inserted "Surface" for this namespace, as else I got this build error, because I simply copied stress_regime.cc:
        // This could for sure be done in a nicer way!
        /**
         * /usr/bin/ld.gold: error: CMakeFiles/aspect.dir/Unity/unity_33.cxx.o: multiple definition of 'aspect::Postprocess::VisualizationPostprocessors::internal::SymmetricTensorImplementation::hybrid(dealii::SymmetricTensor<2, 3, double> const&)'
         * /usr/bin/ld.gold: CMakeFiles/aspect.dir/Unity/unity_34.cxx.o: previous definition here
         */
        namespace SurfaceSymmetricTensorImplementation
        {
          template <int dim>
          void
          tridiagonalize(const dealii::SymmetricTensor<2, dim, double> &A,
                         dealii::Tensor<2, dim, double>                &Q,
                         std::array<double, dim>                       &d,
                         std::array<double, dim - 1> &                  e);

          template <int dim>
          std::array<std::pair<double, Tensor<1, dim, double>>, dim>
                                                            ql_implicit_shifts(const dealii::SymmetricTensor<2, dim, double> &A);
//
//          template <int dim, typename Number>
//          std::array<std::pair<Number, Tensor<1, dim, Number>>, dim>
//            jacobi(dealii::SymmetricTensor<2, dim, Number> A);

          std::array<std::pair<double, Tensor<1, 2, double>>, 2>
                                                          hybrid(const dealii::SymmetricTensor<2, 2, double> &A);

          std::array<std::pair<double, Tensor<1, 3, double>>, 3>
                                                          hybrid(const dealii::SymmetricTensor<2, 3, double> &A);

          std::array<double, 2>
          eigenvalues(const SymmetricTensor<2, 2, double> &T);

          std::array<double, 3>
          eigenvalues(const SymmetricTensor<2, 3, double> &T);

          template <int dim>
          std::array<std::pair<double, Tensor<1, dim, double>>, dim>
                                                            perform_eigenvector_decomposition(const SymmetricTensor<2, dim, double> &T,
                                                                const SurfaceSymmetricTensorEigenvectorMethod method);

          /**
           * A struct that is used to sort arrays of pairs of eign=envalues and
           * eigenvectors. Sorting is performed in descending order of eigenvalue.
           */
          template <int dim>
          struct SortEigenValuesVectors
          {
            using EigValsVecs = std::pair<double, Tensor<1, dim, double>>;
            bool
            operator()(const EigValsVecs &lhs, const EigValsVecs &rhs)
            {
              return lhs.first > rhs.first;
            }
          };
        }
      }
      /**
       * A class that computes a field of horizontal vectors that
       * represent the direction of maximal horizontal compressive
       * stress. For an exact definition, see the documentation of
       * this plugin in the manual.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class SurfaceStressRegime
        : public DataPostprocessor<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>,
          public SurfaceOnlyVisualization<dim>
      {
        public:
          virtual
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double> > &computed_quantities) const;

          /**
           * Return the vector of strings describing the names of the computed
           * quantities. Given the purpose of this class, this is a vector
           * with entries all equal to the name of the plugin.
           */
          virtual std::vector<std::string> get_names () const;

          /**
           * This functions returns information about how the individual
           * components of output files that consist of more than one data set
           * are to be interpreted. The returned value is
           * DataComponentInterpretation::component_is_scalar repeated
           * SymmetricTensor::n_independent_components times. (These
           * components should really be part of a symmetric tensor, but
           * deal.II does not allow marking components as such.)
           */
          virtual
          std::vector<DataComponentInterpretation::DataComponentInterpretation>
          get_data_component_interpretation () const;

          /**
           * Return which data has to be provided to compute the derived
           * quantities. The flags returned here are the ones passed to the
           * constructor of this class.
           */
          virtual UpdateFlags get_needed_update_flags () const;

          /**
           * Return the maximum horizontal compressive stress and the stress regime
           * for a given point.
           */
          std::pair<Tensor<1,dim>, double>
          compute_sigmaH_and_stress_regime(const SymmetricTensor<2,dim> compressive_stress,
                                           const double pressure,
                                           const Tensor<1,dim> vertical_direction,
                                           const std::array<Tensor<1,dim>,dim-1 > orthogonal_directions) const;

        private:
          std::array<std::pair<double, Tensor<1, dim, double> >, dim>
          eigenvectors(const SymmetricTensor<2, dim, double> &,
                       const SurfaceSymmetricTensorEigenvectorMethod) const;


      };
    }
  }
}

#endif
