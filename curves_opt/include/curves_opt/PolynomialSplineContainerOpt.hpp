/*
 * PolynomialSplineContainer.hpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#pragma once

// std
#include <vector>

// curves
#include "curves/PolynomialSplineContainer.hpp"

// numerical optimization
#include <numopt_common/numopt_common.hpp>
#include <numopt_quadprog/ActiveSetFunctionMinimizer.hpp>
#include <numopt_common/ParameterizationIdentity.hpp>
#include <numopt_common/QuadraticProblem.hpp>


namespace curves {

template <int splineOrder_>
class PolynomialSplineContainerOpt : public PolynomialSplineContainer<splineOrder_> {
 public:
  using Base = PolynomialSplineContainer<splineOrder_>;
  using SplineType = typename Base::SplineType;
  using SplineList = std::vector<SplineType>;

  PolynomialSplineContainerOpt();
  virtual ~PolynomialSplineContainerOpt();

  /*!
   * Minimize spline coefficients s.t. position, velocity and acceleration constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the second derivative).
   */
  bool setDataOptimized(
      const std::vector<double>& knotDurations,
      const std::vector<double>& knotPositions,
      double initialVelocity, double initialAcceleration,
      double finalVelocity, double finalAcceleration,
      double weightMinAccel);

  /*!
   * Minimize spline coefficients s.t. position and velocity constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the first derivative).
   */
  bool setDataOptimized(
      const std::vector<double>& knotDurations,
      const std::vector<double>& knotPositions,
      double initialVelocity, double finalVelocity,
      double weightMinAccel);

 protected:
  using MinAccMat = Eigen::Matrix<double, SplineType::coefficientCount-2, SplineType::coefficientCount-2>;

  bool getAccelerationMinimizerBlock(MinAccMat& mat, const double tf) const;

  bool addObjective(
      const std::vector<double>& splineDurations,
      const unsigned int num_coeffs,
      const unsigned int num_splines,
      const double weightMinAccel);

  bool setUpOptimizationMatrices(const unsigned int num_coeffs);

  bool solveQP(Eigen::VectorXd& coeffs);

  // Hessian of quadratic program (Q in x'Qx'l'*x).
  Eigen::MatrixXd hessian_;

  // Linear term of quadratic program (l in x'Qx'l'*x).
  Eigen::VectorXd linearTerm_;

  //! Inequality matrix of quadratic program (A in Ax>=b).
  Eigen::MatrixXd inequalityConstraintJacobian_;

  //! Min values of quadratic program (b in Ax>=b).
  Eigen::VectorXd inequalityConstraintMinValues_;

  //! QP problem.
  std::shared_ptr<numopt_common::QuadraticProblemSolver> minimizer_;
  std::shared_ptr<numopt_common::QuadraticObjectiveFunction> costFunction_;
  std::shared_ptr<numopt_common::LinearFunctionConstraints> functionConstraints_;
  std::shared_ptr<numopt_common::QuadraticProblem> quadraticProblem_;
};

} /* namespace */

#include <curves_opt/PolynomialSplineContainerOpt.tpp>
