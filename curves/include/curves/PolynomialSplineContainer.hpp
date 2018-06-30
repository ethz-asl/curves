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

// std
#include <iostream>
#include <memory>
#include <limits>

// boost
#include <boost/math/special_functions/pow.hpp>

namespace curves {

template <int splineOrder_>
class PolynomialSplineContainer {
 public:
  using SplineType = PolynomialSpline<splineOrder_>;
  using SplineList = std::vector<SplineType>;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  PolynomialSplineContainer();
  virtual ~PolynomialSplineContainer() = default;

  //! Get a pointer to the spline with a given index.
  SplineType* getSpline(int splineIndex);

  //! Get a reference to the spline with a given index.
  const SplineList& getSplines() const;

  //! Update the spline internal time by dt [seconds].
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

  //! Reserve memory for the spline container.
  bool reserveSplines(const unsigned int numSplines);

  //! Clear spline container.
  bool reset();

  //! Set container.
  bool resetTime();

  //! Get total trajectory duration.
  double getContainerDuration() const;

  //! Get currently active time point.
  double getContainerTime() const;

  //! True if splines are empty.
  bool isEmpty() const;

  //! Get the position evaluated at the internal time.
  double getPosition() const;

  //! Get the velocity evaluated at the internal time.
  double getVelocity() const;

  //! Get the acceleration evaluated at the internal time.
  double getAcceleration() const;

  /*! Get the index of the spline active at time t [seconds].
   *  Update timeOffset with the duration of the container at the beginning of the active spline.
   */
  int getActiveSplineIndexAtTime(double t, double& timeOffset) const;

  //! Return spline index at current time validity.
  int getActiveSplineIndex() const;

  //! Get position at time t[seconds];
  double getPositionAtTime(double t) const;

  //! Get velocity at time t[seconds];
  double getVelocityAtTime(double t) const;

  //! Get acceleration at time t[seconds];
  double getAccelerationAtTime(double t) const;

  //! Get position at the end of the spline.
  double getEndPosition() const;

  //! Get velocity at the end of the spline.
  double getEndVelocity() const;

  //! Get acceleration at the end of the spline.
  double getEndAcceleration() const;

  /*!
   * Minimize spline coefficients s.t. position, velocity and acceleration constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the second derivative).
   */
  bool setData(
      const std::vector<double>& knotDurations,
      const std::vector<double>& knotPositions,
      double initialVelocity, double initialAcceleration,
      double finalVelocity, double finalAcceleration);

  /*!
   * Minimize spline coefficients s.t. position and velocity constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the first derivative).
   */
  bool setData(
      const std::vector<double>& knotDurations,
      const std::vector<double>& knotPositions,
      double initialVelocity, double finalVelocity);

  /*!
   * Find linear part of the spline coefficients (a0, a1) s.t. position constraints are satisfied.
   * If the spline order is larger than 1, the remaining spline coefficients are set to zero.
   */
  virtual bool setData(
      const std::vector<double>& knotDurations,
      const std::vector<double>& knotPositions);

  static constexpr double undefinedValue = std::numeric_limits<double>::quiet_NaN();

  bool checkContainer() const;

 protected:
  /*!
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

  //! Equality matrix of quadratic program (A in Ax=b).
  Eigen::MatrixXd equalityConstraintJacobian_;

  //! Equality target values of quatratic program (b in Ax=b).
  Eigen::VectorXd equalityConstraintTargetValues_;
};

} /* namespace */

#include <curves/PolynomialSplineContainer.tpp>
