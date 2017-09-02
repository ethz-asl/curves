/*
 * PolynomialSplineContainer.hpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#pragma once

// Eigen
#include "curves/polynomial_splines.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits>

// numerical optimization
#include <numopt_common/numopt_common.hpp>
#include <numopt_quadprog/ActiveSetFunctionMinimizer.hpp>
#include <numopt_common/ParameterizationIdentity.hpp>
#include "numopt_common/QuadraticProblem.hpp"

// std
#include <iostream>
#include <memory>

// boost
#include <boost/math/special_functions/pow.hpp>

// logger
#include "message_logger/log/log_messages.hpp"

namespace curves {

template <int splineOrder_>
class PolynomialSplineContainer {
 public:

  using SplineType = PolynomialSpline<splineOrder_>;
  using SplineList = std::vector<SplineType>;

  PolynomialSplineContainer();


  virtual ~PolynomialSplineContainer();


  SplineType* getSpline(int splineIndex)
  {
    return &splines_.at(splineIndex);
  }

  const SplineList& getSplines() const {
    return splines_;
  }

  //! move along the trajectory for one sample time.
  bool advance(double dt);

  //! Jump to a specific point in time domain.
  void setContainerTime(double t);

  //! Add a spline to the container.
  template<typename SplineType_>
  bool addSpline(SplineType_&& spline) {
    containerDuration_ += spline.getSplineDuration();
    splines_.emplace_back(std::forward<SplineType_>(spline));
    return true;
  }

  //! Clear spline container
  bool reset();

  //! Set container
  bool resetTime();

  //! Get total trajectory duration.
  double getContainerDuration() const;

  //! Get currently active time point.
  double getContainerTime() const;

  //! True if splines are empty.
  bool isEmpty() const;

  //! Get point in state-space at current time
  double getPosition() const;
  double getVelocity() const;
  double getAcceleration() const;

  //! Get currently active spline index of trajectory.
  int getActiveSplineIndexAtTime(double t, double& timeOffset) const;

  //! Return spline index at current time validity.
  int getActiveSplineIndex() const;

  //! Get pint in state-space at a specific time
  double getPositionAtTime(double t) const;
  double getVelocityAtTime(double t) const;
  double getAccelerationAtTime(double t) const;

  //! Get state-space at final time.
  double getEndPosition() const;
  double getEndVelocity() const;
  double getEndAcceleration() const;

  /*
   * !Minimize spline coefficients s.t. position, velocity and acceleration constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the second derivative).
   * If the spline order is smaller than 5, some boundary constraints are dropped.
   * Notice that the number of junction constraints (and thereby the degree of smoothness at the transition
   * between two adjacent splines) decreases with decreasing spline order.
   * If weightMinAccel is > 0, then the acceleration of the spline segments are minimized. Otherwise, the coefficients
   * spline coefficients are minimized
   */
  bool setData(const std::vector<double>& knotDurations,
               const std::vector<double>& knotPositions,
               double initialVelocity,
               double initialAcceleration,
               double finalVelocity,
               double finalAcceleration,
               double weightMinAccel = -1.0);

  /*
   * !Minimize spline coefficients s.t. position and velocity constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the first derivative).
   * If the spline order is smaller than 3, some boundary constraints are dropped.
   * Notice that the number of junction constraints (and thereby the degree of smoothness at the transition
   * between two adjacent splines) decreases with decreasing spline order.
   * If weightMinAccel is > 0, then the acceleration of the spline segments are minimized. Otherwise, the coefficients
   * spline coefficients are minimized
   */
  bool setData(const std::vector<double>& knotDurations,
               const std::vector<double>& knotPositions,
               double initialVelocity,
               double finalVelocity,
               double weightMinAccel = -1.0);

  /*
   * ! Find linear part of the spline coefficients (a0, a1) s.t. position constraints are satisfied.
   * If the spline order is larger than 1, the remaining spline coefficients are set to zero.
   */
  bool setData(const std::vector<double>& knotDurations,
               const std::vector<double>& knotPositions);

  static constexpr double undefinedValue = std::numeric_limits<double>::quiet_NaN();

 protected:

  using MinAccMat = Eigen::Matrix<double, SplineType::coefficientCount-2, SplineType::coefficientCount-2>;

  /*
   * aijh:
   *  i --> spline id (1,...,n)
   *  j --> spline coefficient aj (a5,...,a1,a0)
   *  h --> dimX, dimY
   *
   * Coefficient vector is:
   *    q = [a15x a14x ... a10x a15y ... a10y a25x ... a20y ... an5x ... an0y]
   */
  inline int getCoeffIndex(const int splineIdx, const int aIdx) const {
    return splineIdx*(splineOrder_+1) + aIdx;
  }

  inline int getSplineColumnIndex(const int splineIdx) const {
    return getCoeffIndex(splineIdx, 0);
  }

  void addInitialConditions(
      const Eigen::VectorXd& initialConditions,
      unsigned int& constraintIdx);

  void addFinalConditions(
      const Eigen::VectorXd& finalConditions,
      unsigned int& constraintIdx,
      const double lastSplineDuration,
      const unsigned int lastSplineId);

  void addJunctionsConditions(
      const std::vector<double>& splineDurations,
      const std::vector<double>& knotPositions,
      unsigned int& constraintIdx,
      const unsigned int num_junctions);

  bool getAccelerationMinimizerBlock(MinAccMat& mat, const double tf) const;

  bool addObjective(
      const std::vector<double>& splineDurations,
      const unsigned int num_coeffs,
      const unsigned int num_splines,
      const double weightMinAccel);

  bool setUpOptimizationMatrices(const unsigned int num_coeffs);

  bool solveQP(Eigen::VectorXd& coeffs);

  bool extractSplineCoefficients(
      const Eigen::VectorXd& coeffs,
      const std::vector<double>& splineDurations,
      const unsigned int num_splines);

  //! Conjunction of smoothly interconnected splines.
  SplineList splines_;

  //! Helper variable.
  double timeOffset_;

  //! Current time of spline conjunction.
  double containerTime_;

  //! Total duration of spline conjunction.
  double containerDuration_;

  //! Spline index currently active.
  int activeSplineIdx_;

  // Hessian of quadratic program (Q in x'Qx'l'*x).
  Eigen::MatrixXd hessian_;

  // Linear term of quadratic program (l in x'Qx'l'*x).
  Eigen::VectorXd linearTerm_;

  //! Equality matrix of quadratic program (A in Ax=b).
  Eigen::MatrixXd equalityConstraintJacobian_;

  //! Equality target values of quatratic program (b in Ax=b).
  Eigen::VectorXd equalityConstraintTargetValues_;

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

#include <curves/PolynomialSplineContainer.tpp>
