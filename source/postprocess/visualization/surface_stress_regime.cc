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


#include <aspect/simulator.h>
#include <aspect/postprocess/visualization/surface_stress_regime.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/utilities.h>

#include <deal.II/physics/transformations.h>

namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      namespace Surfaceinternal
      {
        namespace SurfaceSymmetricTensorImplementation
        {
          template <int dim>
          void
          tridiagonalize(const dealii::SymmetricTensor<2, dim, double> &A,
                         dealii::Tensor<2, dim, double>                &Q,
                         std::array<double, dim>                       &d,
                         std::array<double, dim - 1> &                  e)
          {
            // Create some intermediate storage
            double h, g, omega_inv, K, f;

            // Initialize the transformation matrix as the
            // identity tensor
            Q = dealii::unit_symmetric_tensor<dim, double>();

            // Make the first row and column to be of the
            // desired form
            h = 0.0;
            for (int i = 1; i < dim; i++)
              h += A[0][i] * A[0][i];

            g = 0.0;
            if (A[0][1] > 0.0)
              g = -std::sqrt(h);
            else
              g = std::sqrt(h);
            e[0] = g;

            std::array<double, dim> u;
            for (int i = 1; i < dim; i++)
              {
                u[i] = A[0][i];
                if (i == 1)
                  u[i] -= g;
              }

            std::array<double, dim> q;
            const double omega = h - g * A[0][1];
            if (omega > 0.0)
              {
                omega_inv = 1.0 / omega;
                K         = 0.0;
                for (int i = 1; i < dim; i++)
                  {
                    f = 0.0;
                    for (int j = 1; j < dim; j++)
                      f += A[i][j] * u[j];
                    q[i] = omega_inv * f;
                    K += u[i] * f;
                  }
                K *= 0.5 * omega_inv * omega_inv;

                for (int i = 1; i < dim; i++)
                  q[i] = q[i] - K * u[i];

                d[0] = A[0][0];
                for (int i = 1; i < dim; i++)
                  d[i] = A[i][i] - 2.0 * q[i] * u[i];

                // Store inverse Householder transformation
                // in Q
                for (int j = 1; j < dim; j++)
                  {
                    f = omega_inv * u[j];
                    for (int i = 1; i < dim; i++)
                      Q[i][j] = Q[i][j] - f * u[i];
                  }

                // For dim = 3: Calculate updated A[1][2] and
                // store it in e[1]
                for (int i = 1; i < dim - 1; i++)
                  e[i] = A[i][i + 1] - q[i] * u[i + 1] - u[i] * q[i + 1];
              }
            else
              {
                for (int i = 0; i < dim; i++)
                  d[i] = A[i][i];

                // For dim = 3:
                for (int i = 1; i < dim - 1; i++)
                  e[i] = A[i][i + 1];
              }
          }

          template <int dim>
          std::array<std::pair<double, Tensor<1, dim, double>>, dim>
                                                            ql_implicit_shifts(const dealii::SymmetricTensor<2, dim, double> &A)
          {
            // Transform A to real tridiagonal form by the Householder method:
            // The orthogonal matrix effecting the transformation
            // this will ultimately store the eigenvectors
            dealii::Tensor<2, dim, double> Q;
            // The diagonal elements of the tridiagonal matrix;
            // this will ultimately store the eigenvalues
            std::array<double, dim> w;
            // The off-diagonal elements of the tridiagonal
            std::array<double, dim - 1> ee;
            tridiagonalize<dim>(A, Q, w, ee);

            // Number of iterations
            const unsigned int max_n_it = 30;

            // Transfer the off-diagonal entries to an auxiliary array
            // The third element is used only as temporary workspace
            std::array<double, dim> e;
            for (unsigned int i = 0; i < dim - 1; ++i)
              e[i] = ee[i];

            // Create some intermediate storage
            double g, r, p, f, b, s, c, t;

            // Loop over all off-diagonal elements
            for (int l = 0; l < dim - 1; l++)
              {
                for (unsigned int it = 0; it <= max_n_it; ++it)
                  {
                    // Check for convergence and exit iteration loop
                    // if the off-diagonal element e[l] is zero
                    int m = l;
                    for (; m <= dim - 2; m++)
                      {
                        g = std::abs(w[m]) + std::abs(w[m + 1]);
                        if (std::abs(e[m]) + g == g)
                          break;
                      }
                    if (m == l)
                      break;

                    // Throw if no convergence is achieved within a
                    // stipulated number of iterations
                    if (it == max_n_it)
                      {
                        AssertThrow(
                          false,
                          ExcMessage(
                            "No convergence in iterative QL eigenvector algorithm.")) return std::
                                                                                             array<std::pair<double, Tensor<1, dim, double>>, dim>();
                      }

                    // Calculate the shift..
                    g = (w[l + 1] - w[l]) / (e[l] + e[l]);
                    r = std::sqrt(g * g + 1.0);
                    // .. and then compute g = d_m - k_s for the
                    // plane rotation (Press2007a eq 11.4.22)
                    if (g > 0.0)
                      g = w[m] - w[l] + e[l] / (g + r);
                    else
                      g = w[m] - w[l] + e[l] / (g - r);

                    // Perform plane rotation, as is done in the
                    // standard QL algorithm, followed by Givens
                    // rotations to recover the tridiagonal form
                    s = c = 1.0;
                    p     = 0.0;
                    for (int i = m - 1; i >= l; i--)
                      {
                        f = s * e[i];
                        b = c * e[i];

                        // Branch to recover from underflow
                        if (std::abs(f) > std::abs(g))
                          {
                            c        = g / f;
                            r        = std::sqrt(c * c + 1.0);
                            e[i + 1] = f * r;
                            c *= (s = 1.0 / r);
                          }
                        else
                          {
                            s        = f / g;
                            r        = std::sqrt(s * s + 1.0);
                            e[i + 1] = g * r;
                            s *= (c = 1.0 / r);
                          }

                        g        = w[i + 1] - p;
                        r        = (w[i] - g) * s + 2.0 * c * b;
                        p        = s * r;
                        w[i + 1] = g + p;
                        g        = c * r - b;

                        // Form the eigenvectors
                        for (int k = 0; k < dim; k++)
                          {
                            t           = Q[k][i + 1];
                            Q[k][i + 1] = s * Q[k][i] + c * t;
                            Q[k][i]     = c * Q[k][i] - s * t;
                          }
                      }
                    w[l] -= p;
                    e[l] = g;
                    e[m] = 0.0;
                  }
              }
            // Structure the data to be outputted
            std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs;
            for (unsigned int e = 0; e < dim; ++e)
              {
                eig_vals_vecs[e].first = w[e];

                // The column "e" of Q contains the non-normalized
                // eigenvector associated with the eigenvalue "e"
                for (unsigned int a = 0; a < dim; ++a)
                  {
                    eig_vals_vecs[e].second[a] = Q[a][e];
                  }

                // Normalize
                Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
                eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
              }
            return eig_vals_vecs;
          }

          std::array<std::pair<double, Tensor<1, 2, double>>, 2>
                                                          hybrid(const dealii::SymmetricTensor<2, 2, double> &A)
          {
            const unsigned int dim = 2;

            // Calculate eigenvalues
            const std::array<double, dim> w = eigenvalues(A);

            std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs;

            double t, u; // Intermediate storage
            t = std::abs(w[0]);
            for (unsigned int i = 1; i < dim; ++i)
              {
                u = std::abs(w[i]);
                if (u > t)
                  t = u;
              }

            if (t < 1.0)
              u = t;
            else
              u = t * t;

            // Estimated maximum roundoff error
            const double error =
              256.0 * std::numeric_limits<double>::epsilon() * u * u;

            // Store eigenvalues
            eig_vals_vecs[0].first = w[0];
            eig_vals_vecs[1].first = w[1];

            // Compute eigenvectors
            // http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
            // https://math.stackexchange.com/a/1548616
            if (A[1][0] != 0.0)
              {
                // First eigenvector
                eig_vals_vecs[0].second[0] = w[0] - A[1][1];
                eig_vals_vecs[0].second[1] = A[1][0];

                // Second eigenvector
                eig_vals_vecs[1].second[0] = w[1] - A[1][1];
                eig_vals_vecs[1].second[1] = A[1][0];
              }
            else
              {
                // First eigenvector
                eig_vals_vecs[0].second[0] = w[0];
                eig_vals_vecs[0].second[1] = 0.0;

                // Second eigenvector
                eig_vals_vecs[1].second[0] = 0.0;
                eig_vals_vecs[1].second[1] = w[1];
              }
            // Normalize
            eig_vals_vecs[0].second /= eig_vals_vecs[0].second.norm();
            eig_vals_vecs[1].second /= eig_vals_vecs[1].second.norm();

            // If vectors are nearly linearly dependent, or if there might have
            // been large cancelations in the calculation of A[i][i] - w[0], fall
            // back to QL algorithm
            if (eig_vals_vecs[0].second * eig_vals_vecs[1].second > error)
              {
                return ql_implicit_shifts(A);
              }

            return eig_vals_vecs;
          }


          std::array<std::pair<double, Tensor<1, 3, double>>, 3>
                                                          hybrid(const dealii::SymmetricTensor<2, 3, double> &A)
          {
            const unsigned int dim = 3;
            double norm; // Squared norm or inverse norm of current eigenvector
            double t, u; // Intermediate storage

            // Calculate eigenvalues
            const std::array<double, dim> w = eigenvalues(A);

            t = std::abs(w[0]);
            for (unsigned int i = 1; i < dim; ++i)
              {
                u = std::abs(w[i]);
                if (u > t)
                  t = u;
              }

            if (t < 1.0)
              u = t;
            else
              u = t * t;

            // Estimated maximum roundoff error
            const double error =
              256.0 * std::numeric_limits<double>::epsilon() * u * u;

            // Initialize the transformation matrix as the
            // identity tensor
            dealii::Tensor<2, dim, double> Q;
            Q[0][1] = A[0][1] * A[1][2] - A[0][2] * A[1][1];
            Q[1][1] = A[0][2] * A[0][1] - A[1][2] * A[0][0];
            Q[2][1] = A[0][1] * A[0][1];

            // Calculate first eigenvector by the formula
            //   v[0] = (A - w[0]).e1 x (A - w[0]).e2
            Q[0][0] = Q[0][1] + A[0][2] * w[0];
            Q[1][0] = Q[1][1] + A[1][2] * w[0];
            Q[2][0] = (A[0][0] - w[0]) * (A[1][1] - w[0]) - Q[2][1];
            norm    = Q[0][0] * Q[0][0] + Q[1][0] * Q[1][0] + Q[2][0] * Q[2][0];

            // If vectors are nearly linearly dependent, or if there might have
            // been large cancellations in the calculation of A[i][i] - w[0], fall
            // back to QL algorithm
            // Note that this simultaneously ensures that multiple eigenvalues do
            // not cause problems: If w[0] = w[1], then A - w[0] * I has rank 1,
            // i.e. all columns of A - w[0] * I are linearly dependent.
            if (norm <= error)
              {
                return ql_implicit_shifts(A);
              }
            else // This is the standard branch
              {
                norm = std::sqrt(1.0 / norm);
                for (unsigned j = 0; j < dim; j++)
                  Q[j][0] = Q[j][0] * norm;
              }

            // Calculate second eigenvector by the formula
            //   v[1] = (A - w[1]).e1 x (A - w[1]).e2
            Q[0][1] = Q[0][1] + A[0][2] * w[1];
            Q[1][1] = Q[1][1] + A[1][2] * w[1];
            Q[2][1] = (A[0][0] - w[1]) * (A[1][1] - w[1]) - Q[2][1];
            norm    = Q[0][1] * Q[0][1] + Q[1][1] * Q[1][1] + Q[2][1] * Q[2][1];
            if (norm <= error)
              {
                return ql_implicit_shifts(A);
              }
            else
              {
                norm = std::sqrt(1.0 / norm);
                for (unsigned int j = 0; j < dim; j++)
                  Q[j][1] = Q[j][1] * norm;
              }

            // Calculate third eigenvector according to
            //   v[2] = v[0] x v[1]
            Q[0][2] = Q[1][0] * Q[2][1] - Q[2][0] * Q[1][1];
            Q[1][2] = Q[2][0] * Q[0][1] - Q[0][0] * Q[2][1];
            Q[2][2] = Q[0][0] * Q[1][1] - Q[1][0] * Q[0][1];

            // Structure the data to be outputted
            std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs;
            for (unsigned int e = 0; e < dim; ++e)
              {
                eig_vals_vecs[e].first = w[e];

                // The column "e" of Q contains the non-normalized
                // eigenvector associated with the eigenvalue "e"
                for (unsigned int a = 0; a < dim; ++a)
                  {
                    eig_vals_vecs[e].second[a] = Q[a][e];
                  }

                // Normalize
                Assert(eig_vals_vecs[e].second.norm() != 0.0, ExcDivideByZero());
                eig_vals_vecs[e].second /= eig_vals_vecs[e].second.norm();
              }
            return eig_vals_vecs;
          }


          std::array<double, 2>
          eigenvalues(const SymmetricTensor<2, 2, double> &T)
          {
            const double upp_tri_sq = T[0][1] * T[0][1];
            if (upp_tri_sq == 0.0)
              {
                // The tensor is diagonal
                std::array<double, 2> eig_vals = {{T[0][0], T[1][1]}};

                // Sort from largest to smallest.
                std::sort(eig_vals.begin(), eig_vals.end(), std::greater<double>());
                return eig_vals;
              }
            else
              {
                const double tr_T    = trace(T);
                const double det_T   = determinant(T);
                const double descrim = tr_T * tr_T - 4.0 * det_T;
                Assert(
                  descrim > 0.0,
                  ExcMessage(
                    "The roots of the characteristic polynomial are complex valued."));
                const double sqrt_desc = std::sqrt(descrim);

                const std::array<double, 2> eig_vals =
                {
                  {
                    (0.5 * (tr_T + sqrt_desc)),
                    (0.5 * (tr_T - sqrt_desc))
                  }
                };
                Assert(eig_vals[0] >= eig_vals[1],
                       ExcMessage("The eigenvalue ordering is incorrect."));
                return eig_vals;
              }
          }

          std::array<double, 3>
          eigenvalues(const SymmetricTensor<2, 3, double> &T)
          {
            const double upp_tri_sq =
              T[0][1] * T[0][1] + T[0][2] * T[0][2] + T[1][2] * T[1][2];
            if (upp_tri_sq == 0.0)
              {
                // The tensor is diagonal
                std::array<double, 3> eig_vals = {{T[0][0], T[1][1], T[2][2]}};

                // Sort from largest to smallest.
                std::sort(eig_vals.begin(), eig_vals.end(), std::greater<double>());
                return eig_vals;
              }
            else
              {
                // Perform an affine change to T, and solve a different
                // characteristic equation that has a trigonometric solution.
                // Decompose T = p*B + q*I , and set q = tr(T)/3
                // and p = (tr((T - q.I)^{2})/6)^{1/2} . Then solve the equation
                // 0 = det(\lambda*I - B) = \lambda^{3} - 3*\lambda - det(B)
                // which has the solution
                // \lambda = 2*cos(1/3 * acos(det(B)/2) +2/3*pi*k ); k = 0,1,2
                // when substituting  \lambda = 2.cos(theta) and using trig identities.
                const double tr_T = trace(T);
                const double q    = tr_T / 3.0;
                const double tmp1 = (T[0][0] - q) * (T[0][0] - q) +
                                    (T[1][1] - q) * (T[1][1] - q) +
                                    (T[2][2] - q) * (T[2][2] - q) + 2.0 * upp_tri_sq;
                const double                        p = std::sqrt(tmp1 / 6.0);
                const SymmetricTensor<2, 3, double> B =
                  double(1.0 / p) * (T - q * unit_symmetric_tensor<3, double>());
                const double tmp_2 = determinant(B) / 2.0;

                // The value of tmp_2 should be within [-1,1], however
                // floating point errors might place it slightly outside
                // this range. It is therefore necessary to correct for it.
                // Note: The three results in the conditional may lead to different
                //       number types when using Sacado numbers, so we cast them when
                //       necessary to a consistent result type.
                const double phi =
                  (tmp_2 <= -1.0 ?
                   numbers::PI / 3.0 :
                   (tmp_2 >= 1.0 ?
                    0.0 :
                    std::acos(tmp_2) / 3.0));

                // Due to the trigonometric solution, the computed eigenvalues
                // should be predictably in the order eig1 >= eig2 >= eig3...
                std::array<double, 3> eig_vals =
                {
                  {
                    static_cast<double>(q + 2.0 * p * std::cos(phi)),
                    static_cast<double>(0.0),
                    static_cast<double>(q + 2.0 * p *
                    std::cos(phi + (2.0 / 3.0 * numbers::PI)))
                  }
                };
                // Use the identity tr(T) = eig1 + eig2 + eig3
                eig_vals[1] = tr_T - eig_vals[0] - eig_vals[2];

                // ... however, when equal roots exist then floating point
                // errors may make this no longer be the case.
                // Sort from largest to smallest.
                std::sort(eig_vals.begin(), eig_vals.end(), std::greater<double>());

                return eig_vals;
              }
          }

          template <int dim>
          std::array<std::pair<double, Tensor<1, dim, double>>, dim>
                                                            perform_eigenvector_decomposition(
                                                              const SymmetricTensor<2, dim, double> &T,
                                                              const SurfaceSymmetricTensorEigenvectorMethod method)
          {
            switch (method)
              {
                case SurfaceSymmetricTensorEigenvectorMethod::hybrid:
                  return Surfaceinternal::SurfaceSymmetricTensorImplementation::hybrid(T);
                  break;
                case SurfaceSymmetricTensorEigenvectorMethod::ql_implicit_shifts:
                  //return Surfaceinternal::SurfaceSymmetricTensorImplementation::ql_implicit_shifts(
                  //  T);
                  break;
                case SurfaceSymmetricTensorEigenvectorMethod::jacobi:
                  //  return Surfaceinternal::SurfaceSymmetricTensorImplementation::jacobi(T);
                  break;
                default:
                  break;
              }

            AssertThrow(false, ExcNotImplemented());
            return std::array<std::pair<double, Tensor<1, dim, double>>, dim>();
          }


        } // namespace SurfaceSymmetricTensorImplementation
      } // namespace Surfaceinternal





      template <int dim>
      void
      SurfaceStressRegime<dim>::
      evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                            std::vector<Vector<double> > &computed_quantities) const
      {
        const unsigned int n_quadrature_points = input_data.solution_values.size();
        Assert (computed_quantities.size() == n_quadrature_points,    ExcInternalError());
        // Output is a dim+1 vector of where the first three entries represent the maximum horizontal
        // compressive stress and the fourth has a value of 1, 2 or 3, where 1 is normal faulting, 2 strike-slip and 3 thrust faulting
        Assert ((computed_quantities[0].size() == dim+1),
                ExcInternalError());
        Assert (input_data.solution_values[0].size() == this->introspection().n_components,   ExcInternalError());
        Assert (input_data.solution_gradients[0].size() == this->introspection().n_components,  ExcInternalError());

        MaterialModel::MaterialModelInputs<dim> in(input_data,
                                                   this->introspection());
        MaterialModel::MaterialModelOutputs<dim> out(n_quadrature_points,
                                                     this->n_compositional_fields());

        // Compute the viscosity...
        this->get_material_model().evaluate(in, out);

        // ...and use it to compute the stresses and from that the
        // maximum compressive stress direction
        for (unsigned int q=0; q<n_quadrature_points; ++q)
          {
            const SymmetricTensor<2,dim> strain_rate = in.strain_rate[q];
            const SymmetricTensor<2,dim> compressible_strain_rate
              = (this->get_material_model().is_compressible()
                 ?
                 strain_rate - 1./3 * trace(strain_rate) * unit_symmetric_tensor<dim>()
                 :
                 strain_rate);

            const double eta = out.viscosities[q];

            // first compute the stress tensor, ignoring the pressure
            // for the moment (the pressure has no effect on the
            // direction since it just adds a multiple of the identity
            // matrix to the stress, but because it is large, it may
            // lead to numerical instabilities)
            //
            // note that the *compressive* stress is simply the
            // negative stress
            const SymmetricTensor<2,dim> compressive_stress = -2*eta*compressible_strain_rate;
            const double pressure = in.pressure[q];

            // then find a set of (dim-1) horizontal, unit-length, mutually orthogonal vectors
            const Tensor<1,dim> gravity = this->get_gravity_model().gravity_vector (in.position[q]);
            const Tensor<1,dim> vertical_direction = gravity/gravity.norm();
            std::array<Tensor<1,dim>,dim-1 > orthogonal_directions
              = Utilities::orthogonal_vectors(vertical_direction);
            for (unsigned int i=0; i<orthogonal_directions.size(); ++i)
              orthogonal_directions[i] /= orthogonal_directions[i].norm();

            std::pair<Tensor<1,dim>, double> sigmaH_and_regime = compute_sigmaH_and_stress_regime(compressive_stress, pressure, vertical_direction, orthogonal_directions);

            std::vector<double> sigma_H(dim,0.);

            for (unsigned int i=0; i<dim; ++i)
              {
                sigma_H[i] = (sigmaH_and_regime.first)[i];
                if (std::isnan(sigma_H[i]))
                  {
                    std::cout << "NaN sigma_H component" << std::endl;
                    for (unsigned int j=0; j<dim; ++j)
                      sigma_H[j] = 0.;
                  }
                //computed_quantities[q](i) = maximum_horizontal_compressive_stress[i];
                computed_quantities[q](i) = sigma_H[i];
              }
            computed_quantities[q](dim) = sigmaH_and_regime.second;
          }
      }

      template <int dim>
      std::pair<Tensor<1,dim>, double>
      SurfaceStressRegime<dim>::compute_sigmaH_and_stress_regime(const SymmetricTensor<2,dim> compressive_stress,
                                                                 const double pressure,
                                                                 const Tensor<1,dim> vertical_direction,
                                                                 const std::array<Tensor<1,dim>,dim-1 > orthogonal_directions) const
      {
        double stress_regime = std::numeric_limits<double>::quiet_NaN(), stress_regime_old = std::numeric_limits<double>::quiet_NaN();
        Tensor<1,dim> maximum_horizontal_compressive_stress;
        Tensor<1,dim> sigma_H;
        switch (dim)
          {
            // in 2d, there is only one horizontal direction, and
            // we have already computed it above. give it the
            // length of the compressive_stress (now taking into account the
            // pressure) in this direction
            case 2:
            {
              const double maximum_horizontal_compressive_stress_magnitude = orthogonal_directions[0] *
                                                                             ((compressive_stress
                                                                               -
                                                                               pressure * unit_symmetric_tensor<dim>()) * orthogonal_directions[0]);

              const double vertical_compressive_stress_magnitude = vertical_direction *
                                                                   ((compressive_stress
                                                                     -
                                                                     pressure * unit_symmetric_tensor<dim>()) *
                                                                    vertical_direction);

              // normal faulting
              if (vertical_compressive_stress_magnitude > maximum_horizontal_compressive_stress_magnitude)
                stress_regime = 1.;
              // thrust faulting
              else
                stress_regime = 3.;

              maximum_horizontal_compressive_stress = orthogonal_directions[0] * maximum_horizontal_compressive_stress_magnitude;

              break;
            }

            // in 3d, use the formulas discussed in the
            // documentation of the plugin below
            case 3:
            {
              // Test getting the eigenvectors and -values
              std::array<std::pair<double, Tensor< 1, dim, double> >,dim > eigen_vectors = eigenvectors(compressive_stress, SurfaceSymmetricTensorEigenvectorMethod::hybrid);
              const Tensor<1,dim> S1 = eigen_vectors[0].second;
              const Tensor<1,dim> S2 = eigen_vectors[1].second;
              const Tensor<1,dim> S3 = eigen_vectors[dim-1].second;

              // Compute the plunge and azimuth of the eigenvectors
              // Azimuth is clockwise angle with north, i.e. the angle of
              // the horizontal projection of the eigenvector
              // with a horizontal vector pointing north
              // In cartesian cases, north is a vector pointing in the y-direction.
              // For spherical domains, not implemented.
              // Plunge is the angle with the horizontal plane,
              // e.g. the angle between the vector and its projection on the horizontal
              // plane.
              // 1) project eigenvectors on horizontal plane and normalize
              Tensor<1,dim> S1_hor = S1;
              S1_hor[dim-1] = 0.;
              S1_hor /= S1_hor.norm();
              Tensor<1,dim> S2_hor = S2;
              S2_hor[dim-1] = 0.;
              S2_hor /= S2_hor.norm();
              Tensor<1,dim> S3_hor = S3;
              S3_hor[dim-1] = 0.;
              S3_hor /= S3_hor.norm();

              // 2) compute azimuth through the dot product of a north-pointing
              // unit vector and the projection of the eigenvector. As n.e = ||n|| ||e|| cos(alpha),
              // take the inverse cosine of the dot product divided by the length of the projected
              // eigenvector.
              // Cut-off between -1 and 1 to avoid nans due to round-off error.
              Tensor<1,dim> north;
              north[1] = 1.;
              const double azimuth_S1 = 0.5 * numbers::PI - atan2(S1_hor[1],S1_hor[0]);
              const double azimuth_S2 = 0.5 * numbers::PI - atan2(S2_hor[1],S2_hor[0]);
              const double azimuth_S3 = 0.5 * numbers::PI - atan2(S3_hor[1],S3_hor[0]);

              // 3) compute plunge through the dot product of the projected and true
              // eigenvector. The sign of the plunge is determined by the dot product
              // of the eigenvector and an up-vector, as the plunge is negative for upward
              // pointing eigenvectors.
              Tensor<1,dim> up;
              up[dim-1] = 1.;
              const double plunge_S1 = acos(std::max(-1.,std::min(1., S1 * S1_hor)));
              const double plunge_S2 = acos(std::max(-1.,std::min(1., S2 * S2_hor)));
              const double plunge_S3 = acos(std::max(-1.,std::min(1., S3 * S3_hor)));

              // Distinguish stress regime based on the plunge of the eigenvectors
              // according to the World Stress Map classification as described in
              // Zoback 1992, JGR vol. 97, no. B8.
              // Plunges pl:
              // S1         S2      S3           Regime               sigma_H azimuth
              // >= 52              <=35         normal faulting      S2 azimuth
              // 40<=pl<52          <=20         normal + strikeslip  S3 azimuth + 90 degrees
              // <40        >=45    <=20         strikeslip           S3 azimuth + 90 degrees
              // <=20       >=45    <40          strikeslip           S1 azimuth
              // <=20               40<=pl<52    thrust + strikeslip  S1 azimuth
              // <=35               >=52         thrust faulting      S1 azimuth
              // If the orientation of the eigenvectors does not fit into any of the
              // above regimes, the regime will stay set to NaN.
              const double degree_to_rad = numbers::PI/180.;
              double sigma_H_azimuth = 0.;

              if (plunge_S1 >= 52.*degree_to_rad && plunge_S3 <= 35.*degree_to_rad)
                {
                  stress_regime = 1.;
                  sigma_H_azimuth = azimuth_S2;
                }
              else if (plunge_S1 >= 40.*degree_to_rad &&
                       plunge_S1 <  52.*degree_to_rad &&
                       plunge_S3 <= 35.*degree_to_rad)
                {
                  stress_regime = 2.;
                  sigma_H_azimuth = azimuth_S3 + 90.*degree_to_rad;
                }
              else if (plunge_S1 <  40.*degree_to_rad &&
                       plunge_S2 >= 45.*degree_to_rad &&
                       plunge_S3 <= 20.*degree_to_rad)
                {
                  stress_regime = 3.;
                  sigma_H_azimuth = azimuth_S3 + 90.*degree_to_rad;
                }
              else if (plunge_S1 <= 20.*degree_to_rad &&
                       plunge_S2 >= 45.*degree_to_rad &&
                       plunge_S3 <  40.*degree_to_rad)
                {
                  stress_regime = 4.;
                  sigma_H_azimuth = azimuth_S1;
                }
              else if (plunge_S1 <= 20.*degree_to_rad &&
                       plunge_S3 >= 40.*degree_to_rad &&
                       plunge_S3 <  52.*degree_to_rad)
                {
                  stress_regime = 5.;
                  sigma_H_azimuth = azimuth_S1;
                }
              else if (plunge_S1 <= 35.*degree_to_rad &&
                       plunge_S3 >= 52.*degree_to_rad)
                {
                  stress_regime = 6.;
                  sigma_H_azimuth = azimuth_S1;
                }

              // compute a unit vector in the direction of sigma_H
              Tensor<1,dim> sigma_H_unit;
              sigma_H_unit[0] = std::cos(numbers::PI * 0.5 - sigma_H_azimuth);
              sigma_H_unit[1] = std::sin(numbers::PI * 0.5 - sigma_H_azimuth);

              // compute the maximum horizontal compressive stress
              sigma_H = sigma_H_unit * (sigma_H_unit * ((compressive_stress-pressure * unit_symmetric_tensor<dim>()) * sigma_H_unit));

              const double a = orthogonal_directions[0] *
                               (compressive_stress *
                                orthogonal_directions[0]);
              const double b = orthogonal_directions[1] *
                               (compressive_stress *
                                orthogonal_directions[1]);
              const double c = orthogonal_directions[0] *
                               (compressive_stress *
                                orthogonal_directions[1]);

              // compute the two stationary points of f(alpha)
              const double alpha_1 = 1./2 * std::atan2 (c, a-b);
              const double alpha_2 = alpha_1 + numbers::PI/2;

              // then check the sign of f''(alpha) to determine
              // which of the two stationary points is the maximum
              const double f_double_prime_1 = 2*(b-a)*std::cos(2*alpha_1)
                                              - 2*c*sin(2*alpha_1);
              double alpha;
              if (f_double_prime_1 < 0)
                alpha = alpha_1;
              else
                {
                  Assert (/* f_double_prime_2 = */
                    2*(b-a)*std::cos(2*alpha_2) - 2*c*sin(2*alpha_2) <= 0,
                    ExcInternalError());
                  alpha = alpha_2;
                }

              // then re-assemble the maximum horizontal compressive_stress
              // direction from alpha and the two horizontal
              // directions
              const Tensor<1,dim> n = std::cos(alpha) * orthogonal_directions[0] +
                                      std::sin(alpha) * orthogonal_directions[1];

              // finally compute the actual direction * magnitude,
              // now taking into account the pressure (with the
              // correct sign in front of the pressure for the
              // *compressive* stress)
              //
              // the magnitude is computed as discussed in the
              // description of the plugin below
              const double maximum_horizontal_compressive_stress_magnitude
                = (n * ((compressive_stress
                         -
                         pressure * unit_symmetric_tensor<dim>()) * n));

              const Tensor<1,dim> n_perp = std::sin(alpha) * orthogonal_directions[0] -
                                           std::cos(alpha) * orthogonal_directions[1];

              const double minimum_horizontal_compressive_stress_magnitude
                = (n_perp * ((compressive_stress
                              -
                              pressure * unit_symmetric_tensor<dim>()) * n_perp));

              const double vertical_compressive_stress_magnitude
                = (vertical_direction * ((compressive_stress
                                          -
                                          pressure * unit_symmetric_tensor<dim>()) * vertical_direction));

              // normal faulting
              if (vertical_compressive_stress_magnitude > maximum_horizontal_compressive_stress_magnitude &&
                  maximum_horizontal_compressive_stress_magnitude > minimum_horizontal_compressive_stress_magnitude)
                stress_regime_old = 1.;
              // strike-slip faulting
              else if (maximum_horizontal_compressive_stress_magnitude > vertical_compressive_stress_magnitude &&
                       vertical_compressive_stress_magnitude > minimum_horizontal_compressive_stress_magnitude)
                stress_regime_old = 3.;
              // thrust faulting
              else if (maximum_horizontal_compressive_stress_magnitude > minimum_horizontal_compressive_stress_magnitude &&
                       minimum_horizontal_compressive_stress_magnitude > vertical_compressive_stress_magnitude)
                stress_regime_old = 6.;

              maximum_horizontal_compressive_stress = n * (maximum_horizontal_compressive_stress_magnitude - minimum_horizontal_compressive_stress_magnitude);

              break;
            }


            default:
              Assert (false, ExcNotImplemented());
          }

        return std::make_pair(sigma_H, stress_regime);

      }

      template <int dim>
      //std::array<std::pair< Number, Tensor<1, dim, Number> >,std::integral_constant<int, dim>::value>
      std::array<std::pair<double, Tensor<1, dim, double>>,dim>
                                                        SurfaceStressRegime<dim>::eigenvectors(const SymmetricTensor<2, dim, double> &T,
                                                            const SurfaceSymmetricTensorEigenvectorMethod method) const
      {
        std::array<std::pair<double, Tensor<1, dim, double>>, dim> eig_vals_vecs = Surfaceinternal::SurfaceSymmetricTensorImplementation::
                                                                                   perform_eigenvector_decomposition(T, method);

        // Sort in descending order before output.
        std::sort(
          eig_vals_vecs.begin(),
          eig_vals_vecs.end(),
          Surfaceinternal::SurfaceSymmetricTensorImplementation::SortEigenValuesVectors<dim>());
        return eig_vals_vecs;
      }

      template <int dim>
      std::vector<std::string>
      SurfaceStressRegime<dim>::get_names () const
      {
        std::vector<std::string> names(dim+1);
        names[0] = "surface sigma_H_x";
        names[dim-1] = "surface sigma_H_z";
        names[1] = "surface sigma_H_y";
        names[dim] = "surface stress_regime";
        return names;
      }


      template <int dim>
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      SurfaceStressRegime<dim>::get_data_component_interpretation () const
      {
        std::vector<DataComponentInterpretation::DataComponentInterpretation> interpretation
        (dim+1,DataComponentInterpretation::component_is_part_of_vector);
        interpretation[dim] = DataComponentInterpretation::component_is_scalar;
        return interpretation;
      }



      template <int dim>
      UpdateFlags
      SurfaceStressRegime<dim>::get_needed_update_flags () const
      {
        return update_gradients | update_values; // | update_q_points;
      }

    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      ASPECT_REGISTER_VISUALIZATION_POSTPROCESSOR(SurfaceStressRegime,
                                                  "surface stress regime",
                                                  "A plugin that computes the direction of the "
                                                  "maximum horizontal component of the compressive stress as a vector "
                                                  "field, scaled with a value that indicates the principle mode of deformation. "
                                                  "A value of 1 indicates normal faulting, 2 strike-slip and 3 thrust faulting. "
                                                  "Recall that the "
                                                  "\\textit{compressive} stress is simply the negative stress, "
                                                  "$\\sigma_c=-\\sigma=-\\left["
                                                  "     2\\eta (\\varepsilon(\\mathbf u)"
                                                  "             - \\frac 13 (\\nabla \\cdot \\mathbf u) I)"
                                                  "     + pI\\right]$."
                                                  "\n\n"
                                                  "Following \\cite{LundTownend07}, we define the maximum horizontal "
                                                  "stress direction as that \\textit{horizontal} direction "
                                                  "$\\mathbf n$ that maximizes $\\mathbf n^T \\sigma_c \\mathbf n$. We "
                                                  "call a vector \\textit{horizontal} if it is perpendicular to the "
                                                  "gravity vector $\\mathbf g$."
                                                  "\n\n"
                                                  "In two space dimensions, $\\mathbf n$ is simply a vector that "
                                                  "is horizontal (we choose one of the two possible choices). "
                                                  "This direction is then scaled by the size of the horizontal stress "
                                                  "in this direction, i.e., the plugin outputs the vector "
                                                  "$\\mathbf w = (\\mathbf n^T \\sigma_c \\mathbf n) \\; \\mathbf n$."
                                                  "\n\n"
                                                  "In three space dimensions, given two horizontal, perpendicular, "
                                                  "unit length, but otherwise arbitrarily chosen vectors "
                                                  "$\\mathbf u,\\mathbf v$, we can express "
                                                  "$\\mathbf n = (\\cos \\alpha)\\mathbf u + (\\sin\\alpha)\\mathbf v$ "
                                                  "where $\\alpha$ maximizes the expression "
                                                  "\\begin{align*}"
                                                  "  f(\\alpha) = \\mathbf n^T \\sigma_c \\mathbf n"
                                                  "  = (\\mathbf u^T \\sigma_c \\mathbf u)(\\cos\\alpha)^2"
                                                  "    +2(\\mathbf u^T \\sigma_c \\mathbf v)(\\cos\\alpha)(\\sin\\alpha)"
                                                  "    +(\\mathbf v^T \\sigma_c \\mathbf v)(\\sin\\alpha)^2."
                                                  "\\end{align*}"
                                                  "\n\n"
                                                  "The maximum of $f(\\alpha)$ is attained where $f'(\\alpha)=0$. "
                                                  "Evaluating the derivative and using trigonometric identities, "
                                                  "one finds that $\\alpha$ has to satisfy the equation "
                                                  "\\begin{align*}"
                                                  "  \\tan(2\\alpha) = \\frac{\\mathbf u^T \\sigma_c \\mathbf v}"
                                                  "                          {\\mathbf u^T \\sigma_c \\mathbf u "
                                                  "                           - \\mathbf v^T \\sigma_c \\mathbf v}."
                                                  "\\end{align*}"
                                                  "Since the transform $\\alpha\\mapsto\\alpha+\\pi$ flips the "
                                                  "direction of $\\mathbf n$, we only need to seek a solution "
                                                  "to this equation in the interval $\\alpha\\in[0,\\pi)$. "
                                                  "These are given by "
                                                  "$\\alpha_1=\\frac 12 \\arctan \\frac{\\mathbf u^T \\sigma_c "
                                                  "\\mathbf v}{\\mathbf u^T \\sigma_c \\mathbf u - "
                                                  "\\mathbf v^T \\sigma_c \\mathbf v}$ and "
                                                  "$\\alpha_2=\\alpha_1+\\frac{\\pi}{2}$, one of which will "
                                                  "correspond to a minimum and the other to a maximum of "
                                                  "$f(\\alpha)$. One checks the sign of "
                                                  "$f''(\\alpha)=-2(\\mathbf u^T \\sigma_c \\mathbf u - "
                                                  "\\mathbf v^T \\sigma_c \\mathbf v)\\cos(2\\alpha) "
                                                  "- 2 (\\mathbf u^T \\sigma_c \\mathbf v) \\sin(2\\alpha)$ for "
                                                  "each of these to determine the $\\alpha$ that maximizes "
                                                  "$f(\\alpha)$, and from this immediately arrives at the correct "
                                                  "form for the maximum horizontal stress $\\mathbf n$."
                                                  "\n\n"
                                                  "The stress regime is determined based on the magnitudes of the "
                                                  "maximum ($\\sigma_H$) and minimum ($\\sigma_h$) horizontal stress "
                                                  "and the vertical stress, "
                                                  "where $\\sigma_v>\\sigma_H>\\sigma_h$ defines normal faulting, "
                                                  "$\\sigma_H>\\sigma_v>\\sigma_h$ strike-slip faulting and "
                                                  "$\\sigma_H>\\sigma_h>\\sigma_v$ thrust faulting. "
                                                 )
    }
  }
}
