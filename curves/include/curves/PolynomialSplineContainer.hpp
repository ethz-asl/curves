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


  ~PolynomialSplineContainer();


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

  //! Add (push back) a spline.
  bool addSpline(const SplineType& spline);

  //! Add (push back) a spline.
  bool addSpline(SplineType&& spline);

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


  /*
   * aijh:
   *  i --> spline id (1,...,n)
   *  j --> spline coefficient aj (a5,...,a1,a0)
   *  h --> dimX, dimY
   *
   * Coefficient vector is:
   *    q = [a15x a14x ... a10x a15y ... a10y a25x ... a20y ... an5x ... an0y]
   */
  int getCoeffIndex(int splineIdx, int aIdx) const;

  int getSplineColumnIndex(int splineIdx) const;

  void addInitialConditions(const Eigen::VectorXd& initialConditions,
                            unsigned int& constraintIdx);

  void addFinalConditions(const Eigen::VectorXd& finalConditions,
                          unsigned int& constraintIdx,
                          double lastSplineDuration,
                          unsigned int lastSplineId);

  void addJunctionsConditions(const std::vector<double>& splineDurations,
                              const std::vector<double>& knotPositions,
                              unsigned int& constraintIdx,
                              unsigned int num_junctions);

  bool getAccelerationMinimizerBlock(Eigen::MatrixXd& mat, double tf) const;

  bool addObjective(
      const std::vector<double>& splineDurations,
      unsigned int num_coeffs,
      unsigned int num_splines,
      double weightMinAccel);

  bool setUpOptimizationMatrices(unsigned int num_coeffs);

  bool solveQP(Eigen::VectorXd& coeffs, unsigned int num_coeffs);

  bool extractSplineCoefficients(const Eigen::VectorXd& coeffs,
                                 const std::vector<double>& splineDurations,
                                 unsigned int num_splines,
                                 unsigned int num_coeffs_spline);

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

  //! Min values of quatratic program (b in Ax>=b).
  Eigen::VectorXd inequalityConstraintMinValues_;

  //! QP problem.
  std::shared_ptr<numopt_common::QuadraticProblemSolver> minimizer_;
  std::shared_ptr<numopt_common::QuadraticObjectiveFunction> costFunction_;
  std::shared_ptr<numopt_common::LinearFunctionConstraints> functionConstraints_;
  std::shared_ptr<numopt_common::QuadraticProblem> quadraticProblem_;
};


template <int splineOrder_>
PolynomialSplineContainer<splineOrder_>::PolynomialSplineContainer():
    timeOffset_(0.0),
    containerTime_(0.0),
    containerDuration_(0.0),
    activeSplineIdx_(0),
    hessian_(),
    linearTerm_(),
    equalityConstraintJacobian_(),
    equalityConstraintTargetValues_(),
    inequalityConstraintJacobian_(),
    inequalityConstraintMinValues_()
{
  // Make sure that the container is correctly emptied
  reset();

  minimizer_.reset(new numopt_quadprog::ActiveSetFunctionMinimizer());
  costFunction_.reset(new numopt_common::QuadraticObjectiveFunction());
  functionConstraints_.reset(new numopt_common::LinearFunctionConstraints());
  quadraticProblem_.reset(new numopt_common::QuadraticProblem(costFunction_, functionConstraints_));
}

template <int splineOrder_>
PolynomialSplineContainer<splineOrder_>::~PolynomialSplineContainer()
{
  // TODO Auto-generated destructor stub
}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::advance(double dt)
{
  if (splines_.empty() || containerTime_ >= containerDuration_ || activeSplineIdx_ == splines_.size()) {
    return false;
  }

  // advance in time.
  containerTime_ += dt;

  // check if spline index needs to be increased.
  if ((containerTime_ - timeOffset_ >= splines_[activeSplineIdx_].getSplineDuration())) {
    if (activeSplineIdx_ < (splines_.size() - 1)) {
      timeOffset_ += splines_[activeSplineIdx_].getSplineDuration();
    }
    activeSplineIdx_++;
  }

  return true;
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::setContainerTime(double t)
{
  containerTime_ = t;
  double timeOffset;
  activeSplineIdx_ = getActiveSplineIndexAtTime(t, timeOffset);
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::addSpline(const SplineType& spline)
{
  splines_.push_back(spline);
  containerDuration_ += spline.getSplineDuration();
  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::addSpline(SplineType&& spline) {
  containerDuration_ += spline.getSplineDuration();
  splines_.emplace_back(spline);
  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::reset()
{
  splines_.clear();
  activeSplineIdx_ = 0;
  containerDuration_ = 0.0;
  resetTime();
  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::resetTime()
{
  timeOffset_ = 0.0;
  containerTime_ = 0.0;
  activeSplineIdx_ = 0;
  return true;
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getContainerDuration() const
{
  return containerDuration_;
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getContainerTime() const
{
  return containerTime_;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::isEmpty() const
{
  return splines_.empty();
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getPosition() const
{
  if (splines_.empty()) return 0.0;
  if (activeSplineIdx_ == splines_.size())
    return splines_.at(activeSplineIdx_ - 1).getPositionAtTime(containerTime_ - timeOffset_);
  return splines_.at(activeSplineIdx_).getPositionAtTime(containerTime_ - timeOffset_);
}


template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getVelocity() const {
  if (splines_.empty()) return 0.0;
  if (activeSplineIdx_ == splines_.size())
    return splines_.at(activeSplineIdx_ - 1).getVelocityAtTime(containerTime_ - timeOffset_);
  return splines_.at(activeSplineIdx_).getVelocityAtTime(containerTime_ - timeOffset_);
}


template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getAcceleration() const
{
  if (splines_.empty()) return 0.0;
  if (activeSplineIdx_ == splines_.size())
    return splines_.at(activeSplineIdx_ - 1).getAccelerationAtTime(containerTime_ - timeOffset_);
  return splines_.at(activeSplineIdx_).getAccelerationAtTime(containerTime_ - timeOffset_);
}

template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getActiveSplineIndexAtTime(double t, double& timeOffset) const {
  if (splines_.empty()) return -1;
  timeOffset = 0.0;

  for (size_t i = 0; i < splines_.size(); i++) {
    if ((t - timeOffset < splines_[i].getSplineDuration())) {
      return i;
    }
    if (i < (splines_.size() - 1)) {
      timeOffset += splines_[i].getSplineDuration();
    }
  }

  return (splines_.size() - 1);
}

template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getActiveSplineIndex() const
{
  return activeSplineIdx_;
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getPositionAtTime(double t) const
{
  double timeOffset = 0.0;
  int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);

  if (activeSplineIdx < 0) {
    return splines_.at(0).getPositionAtTime(0.0);
  }
  if (activeSplineIdx == splines_.size()) {
    return splines_.at(activeSplineIdx - 1).getPositionAtTime(
        splines_.at(activeSplineIdx - 1).getSplineDuration());
  }
  return splines_.at(activeSplineIdx).getPositionAtTime(t - timeOffset);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getVelocityAtTime(double t) const
{
  double timeOffset = 0.0;
  int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);
  if (activeSplineIdx < 0) {
    return splines_.at(0).getVelocityAtTime(0.0);
  }
  if (activeSplineIdx == splines_.size()) {
    return splines_.at(activeSplineIdx - 1).getVelocityAtTime(
        splines_.at(activeSplineIdx - 1).getSplineDuration());
  }
  return splines_.at(activeSplineIdx).getVelocityAtTime(t - timeOffset);
}


template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getAccelerationAtTime(double t) const
{
  double timeOffset = 0.0;
  int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);
  if (activeSplineIdx < 0) {
    return splines_.at(0).getAccelerationAtTime(0.0);
  }
  if (activeSplineIdx == splines_.size()) {
    return splines_.at(activeSplineIdx - 1).getAccelerationAtTime(
        splines_.at(activeSplineIdx - 1).getSplineDuration());
  }
  return splines_.at(activeSplineIdx).getAccelerationAtTime(t - timeOffset);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getEndPosition() const
{
  double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
  return splines_.at(splines_.size() - 1).getPositionAtTime(lastSplineDuration);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getEndVelocity() const
{
  double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
  return splines_.at(splines_.size() - 1).getVelocityAtTime(lastSplineDuration);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getEndAcceleration() const
{
  double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
  return splines_.at(splines_.size() - 1).getAccelerationAtTime(lastSplineDuration);
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setData(const std::vector<double>& knotDurations,
             const std::vector<double>& knotPositions,
             double initialVelocity,
             double initialAcceleration,
             double finalVelocity,
             double finalAcceleration,
             double weightMinAccel) {

  bool success = true;

  success &= reset();

  // Set up optimization parameters
  const unsigned int num_splines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  const unsigned int num_coeffs = num_splines*num_coeffs_spline;
  const unsigned int num_junctions = num_splines-1;

  if (num_splines<1) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Not sufficient knot points are available!");
    return false;
  }

  // total number of constraints
  constexpr unsigned int num_initial_constraints = 3;   // pos, vel, accel
  constexpr unsigned int num_final_constraints = 3;     // pos, vel, accel
  constexpr unsigned int num_constraint_junction = 4;   // pos (2x), vel, accel
  const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
  const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

  // drop constraints if necessary
  if (num_constraints>num_coeffs) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Number of equality constraints is larger than number of coefficients. Drop acceleration constraints!");
    return setData(knotDurations, knotPositions, initialVelocity, finalVelocity, weightMinAccel);
  }

  // Vector containing durations of splines
  std::vector<double> splineDurations(num_splines);
  for (unsigned int splineId=0; splineId<num_splines; splineId++) {
    splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

    if (splineDurations[splineId]<=0.0) {
      MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Invalid spline duration!");
      return false;
    }
  }

  // Initialize Equality matrices
  equalityConstraintJacobian_.setZero(num_constraints, num_coeffs);
  equalityConstraintTargetValues_.setZero(num_constraints);
  unsigned int constraintIdx = 0;

  // Initial conditions
  Eigen::VectorXd initialConditions(num_initial_constraints);
  initialConditions(0) = knotPositions.front();
  initialConditions(1) = initialVelocity;
  initialConditions(2) = initialAcceleration;
  addInitialConditions(initialConditions, constraintIdx);

  // Final conditions
  Eigen::VectorXd finalConditions(num_final_constraints);
  finalConditions(0) = knotPositions.back();
  finalConditions(1) = finalVelocity;
  finalConditions(2) = finalAcceleration;
  addFinalConditions(finalConditions, constraintIdx,  splineDurations.back(), num_junctions);

  // Junction conditions
  addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

  if (num_constraints != constraintIdx) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Wrong number of equality constraints!");
    return false;
  }

  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(num_coeffs);

  // Minimize spline coefficients
  if (weightMinAccel<=0.0 || num_constraints==num_coeffs) {
    coeffs = equalityConstraintJacobian_.colPivHouseholderQr().solve(equalityConstraintTargetValues_);
  }

  // Minimize acceleration
  else {
    success &= addObjective(splineDurations, num_coeffs, num_splines, weightMinAccel);
    success &= setUpOptimizationMatrices(num_coeffs);
    success &= solveQP(coeffs, num_coeffs);
  }

  // Extract spline coefficients and add splines
  success &= extractSplineCoefficients(coeffs, splineDurations, num_splines, num_coeffs_spline);

  return success;

}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setData(const std::vector<double>& knotDurations,
             const std::vector<double>& knotPositions,
             double initialVelocity,
             double finalVelocity,
             double weightMinAccel) {

  bool success = true;

  success &= reset();

  // Set up optimization parameters
  const unsigned int num_splines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  const unsigned int num_coeffs = num_splines*num_coeffs_spline;
  const unsigned int num_junctions = num_splines-1;

  if (num_splines<1) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Not sufficient knot points are available!");
    return false;
  }

  // total number of constraints
  constexpr unsigned int num_initial_constraints = 2;   // pos, vel
  constexpr unsigned int num_final_constraints = 2;     // pos, vel
  constexpr unsigned int num_constraint_junction = 3;   // pos (2x), vel
  const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
  const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

  // drop constraints if necessary
  if (num_constraints>num_coeffs) {
    std::cout << "[PolynomialSplineContainer::setData] Number of equality constraints is larger than number of coefficients. Drop velocity constraints!\n";
    return setData(knotDurations, knotPositions);
  }

  // Vector containing durations of splines
  std::vector<double> splineDurations(num_splines);
  for (unsigned int splineId=0; splineId<num_splines; splineId++) {
    splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

    if (splineDurations[splineId]<=0.0) {
      std::cout << "[PolynomialSplineContainer::setData] Invalid spline duration!\n";
      return false;
    }
  }

  // Initialize Equality matrices
  equalityConstraintJacobian_.setZero(num_constraints, num_coeffs);
  equalityConstraintTargetValues_.setZero(num_constraints);
  unsigned int constraintIdx = 0;

  // Initial conditions
  Eigen::VectorXd initialConditions(num_initial_constraints);
  initialConditions(0) = knotPositions.front();
  initialConditions(1) = initialVelocity;
  addInitialConditions(initialConditions, constraintIdx);

  // Final conditions
  Eigen::VectorXd finalConditions(num_final_constraints);
  finalConditions(0) = knotPositions.back();
  finalConditions(1) = finalVelocity;
  addFinalConditions(finalConditions, constraintIdx,  splineDurations.back(), num_junctions);

  // Junction conditions
  addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

  if (num_constraints!=constraintIdx) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Wrong number of equality constraints!");
    return false;
  }

  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(num_coeffs);

  // Minimize spline coefficients
  if (weightMinAccel<=0.0 || num_constraints==num_coeffs) {
    coeffs = equalityConstraintJacobian_.colPivHouseholderQr().solve(equalityConstraintTargetValues_);
  }

  // Minimize acceleration
  else {
    success &= addObjective(splineDurations, num_coeffs, num_splines, weightMinAccel);
    success &= setUpOptimizationMatrices(num_coeffs);
    success &= solveQP(coeffs, num_coeffs);
  }

  // Extract spline coefficients and add splines
  success &= extractSplineCoefficients(coeffs, splineDurations, num_splines, num_coeffs_spline);

  return success;

}



template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setData(const std::vector<double>& knotDurations,
             const std::vector<double>& knotPositions) {

  bool success = true;
  const int num_splines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  typename SplineType::SplineCoefficients coefficients;

  if (splineOrder_==0 || num_splines<1) {
    return false;
  }

  std::fill(coefficients.begin(), coefficients.end(), 0.0);

  success &= reset();

  for (unsigned int splineId = 0; splineId<num_splines; splineId++) {
    const double duration = knotDurations[splineId+1]-knotDurations[splineId];

    if (duration<=0.0) {
      return false;
    }

    coefficients[num_coeffs_spline-1] = knotPositions[splineId]; //a0
    coefficients[num_coeffs_spline-2] = (knotPositions[splineId+1]-knotPositions[splineId]) / duration; // a1
    success &= this->addSpline(SplineType(coefficients, duration));
  }

  return success;

}


template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getCoeffIndex(int splineIdx, int aIdx) const {
  return splineIdx*(splineOrder_+1) + aIdx;
}

template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getSplineColumnIndex(int splineIdx) const
{
  return getCoeffIndex(splineIdx, 0);
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::addInitialConditions(const Eigen::VectorXd& initialConditions,
                          unsigned int& constraintIdx) {

  // time container
  typename SplineType::EigenTimeVectorType timeVec;

  // initial position
  if (initialConditions.size()>0) {
    SplineType::getTimeVector(timeVec, 0.0);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = initialConditions(0);
    constraintIdx++;
  }

  // initial velocity
  if (initialConditions.size()>1) {
    SplineType::getdTimeVector(timeVec, 0.0);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = initialConditions(1);
    constraintIdx++;
  }

  // initial acceleration
  if (initialConditions.size()>2) {
    SplineType::getddTimeVector(timeVec, 0.0);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = initialConditions(2);
    constraintIdx++;
  }
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::addFinalConditions(const Eigen::VectorXd& finalConditions,
                        unsigned int& constraintIdx,
                        double lastSplineDuration,
                        unsigned int lastSplineId) {

  // time container
  typename SplineType::EigenTimeVectorType timeVec0;


  // initial position
  if (finalConditions.size()>0) {
    SplineType::getTimeVector(timeVec0, lastSplineDuration);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec0;
    equalityConstraintTargetValues_(constraintIdx) = finalConditions(0);
    constraintIdx++;
  }

  // initial velocity
  if (finalConditions.size()>1) {
    SplineType::getdTimeVector(timeVec0, lastSplineDuration);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec0;
    equalityConstraintTargetValues_(constraintIdx) = finalConditions(1);
    constraintIdx++;
  }

  // initial acceleration
  if (finalConditions.size()>2) {
    SplineType::getddTimeVector(timeVec0, lastSplineDuration);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec0;
    equalityConstraintTargetValues_(constraintIdx) = finalConditions(2);
    constraintIdx++;
  }
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::addJunctionsConditions(const std::vector<double>& splineDurations,
                            const std::vector<double>& knotPositions,
                            unsigned int& constraintIdx,
                            unsigned int num_junctions) {

  // time containers
  typename SplineType::EigenTimeVectorType timeVec0, dTimeVec0, ddTimeVec0;
  typename SplineType::EigenTimeVectorType timeVecTf, dTimeVecTf, ddTimeVecTf;

  // get time container at zero time
  SplineType::getTimeVector(timeVec0, 0.0);
  SplineType::getdTimeVector(dTimeVec0, 0.0);
  SplineType::getddTimeVector(ddTimeVec0, 0.0);


  for (unsigned int splineId=0; splineId<num_junctions; splineId++) {
    const unsigned int nextSplineId = splineId+1;

    // get time container at spline duration
    SplineType::getTimeVector(timeVecTf, splineDurations[splineId]);
    SplineType::getdTimeVector(dTimeVecTf, splineDurations[splineId]);
    SplineType::getddTimeVector(ddTimeVecTf, splineDurations[splineId]);

    // smooth position transition with fixed positions
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  timeVecTf;
    equalityConstraintTargetValues_(constraintIdx) = knotPositions[nextSplineId];
    constraintIdx++;

    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) =  timeVec0;
    equalityConstraintTargetValues_(constraintIdx) = knotPositions[nextSplineId];
    constraintIdx++;

    // smooth velocity transition
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  dTimeVecTf;
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) = -dTimeVec0;
    equalityConstraintTargetValues_(constraintIdx) = 0.0;
    constraintIdx++;

    // smooth acceleration transition
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  ddTimeVecTf;
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) = -ddTimeVec0;
    equalityConstraintTargetValues_(constraintIdx) = 0.0;
    constraintIdx++;
  }
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::getAccelerationMinimizerBlock(Eigen::MatrixXd& mat, double tf) const {
  mat.resize(splineOrder_-1, splineOrder_-1);

  if (splineOrder_ == 5) {
    const double tf2 = boost::math::pow<2>(tf);
    const double tf3 = tf2*tf;
    const double tf4 = tf3*tf;
    const double tf5 = tf4*tf;
    const double tf6 = tf5*tf;
    const double tf7 = tf6*tf;

    mat << 400.0/7.0*tf7, 40.0*tf6,       24.0*tf5,       10.0*tf4,
           40.0*tf6,      28.8*tf5,       18.0*tf4,       8.0*tf3,
           24.0*tf5,      18.0*tf4,       12.0*tf3,       6.0*tf2,
           10.0*tf4,      8.0*tf3,        6.0*tf2,        4.0*tf;
  }

  else if (splineOrder_ == 4) {
    const double tf2 = boost::math::pow<2>(tf);
    const double tf3 = tf2*tf;
    const double tf4 = tf3*tf;
    const double tf5 = tf4*tf;

    mat << 28.8*tf5,       18.0*tf4,       8.0*tf3,
           18.0*tf4,       12.0*tf3,       6.0*tf2,
           8.0*tf3,        6.0*tf2,        4.0*tf;

  }

  else if (splineOrder_ == 3) {
    const double tf2 = boost::math::pow<2>(tf);
    const double tf3 = tf2*tf;

    mat << 12.0*tf3,       6.0*tf2,
           6.0*tf2,        4.0*tf;
  }

  else if (splineOrder_ == 2) {
    mat << 4.0*tf;
  }

  else {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData::getAccelerationMinimizerBlock] Function has not been implemented so far.");
    return false;
  }

  return true;
}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::addObjective(
    const std::vector<double>& splineDurations,
    unsigned int num_coeffs,
    unsigned int num_splines,
    double weightMinAccel) {

  bool success = true;

  // Objective -> minimize acceleration along trajectory
  hessian_ = Eigen::MatrixXd::Zero(num_coeffs,num_coeffs);
  linearTerm_ = Eigen::VectorXd::Zero(num_coeffs);

  for (unsigned int splineId=0; splineId<num_splines; splineId++) {

    // Minimize acceleration
    Eigen::MatrixXd coreMatrix;
    success &= getAccelerationMinimizerBlock(coreMatrix, splineDurations[splineId]);
    hessian_.block(getSplineColumnIndex(splineId),  getSplineColumnIndex(splineId), splineOrder_-1, splineOrder_-1) = coreMatrix*weightMinAccel;

    // Regularization for position and velocity -> cfMat must be positive definite.
    const unsigned int accelIdx = getCoeffIndex(splineId, splineOrder_-1);
    hessian_.block(accelIdx, accelIdx, 2, 2).noalias() += 1e-7*Eigen::Matrix2d::Identity();
  }

  return success;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setUpOptimizationMatrices(unsigned int num_coeffs) {

  // Inequality constraints -> We don't use them
  inequalityConstraintJacobian_.setZero(0, num_coeffs);
  inequalityConstraintMinValues_.setZero(0);

  // Add quadratic problem
  costFunction_->setGlobalHessian(hessian_.sparseView());
  costFunction_->setLinearTerm(linearTerm_);
  functionConstraints_->setGlobalEqualityConstraintJacobian(equalityConstraintJacobian_.sparseView());
  functionConstraints_->setEqualityConstraintTargetValues(equalityConstraintTargetValues_);
  functionConstraints_->setGlobalInequalityConstraintJacobian(inequalityConstraintJacobian_.sparseView());
  functionConstraints_->setInequalityConstraintMinValues(inequalityConstraintMinValues_);
  functionConstraints_->setInequalityConstraintMaxValues(inequalityConstraintMinValues_);


  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::solveQP(Eigen::VectorXd& coeffs, unsigned int num_coeffs) {
  bool success = true;

  double cost = 0.0;
  numopt_common::ParameterizationIdentity params(num_coeffs);
  params.getParams() = coeffs;
  success &= minimizer_->minimize(quadraticProblem_.get(), params, cost);

  if (!success) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData::solveQP] Failed to solve optimization!");
  } else {
    coeffs = params.getParams();
  }

  return success;
}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::extractSplineCoefficients(
                               const Eigen::VectorXd& coeffs,
                               const std::vector<double>& splineDurations,
                               unsigned int num_splines,
                               unsigned int num_coeffs_spline) {
  SplineType spline;
  typename SplineType::SplineCoefficients coefficients;

  for (unsigned int splineId = 0; splineId<num_splines; splineId++) {
    Eigen::Map<Eigen::VectorXd>(coefficients.data(), num_coeffs_spline, 1) = coeffs.segment(getSplineColumnIndex(splineId), num_coeffs_spline);
    spline.setCoefficientsAndDuration(coefficients, splineDurations[splineId]);
    this->addSpline(spline);
  }

  return true;
}




} /* namespace */
