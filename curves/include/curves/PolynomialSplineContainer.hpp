/*
 * PolynomialSplineContainer.hpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#pragma once

#include "curves/polynomial_splines.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits>

namespace curves {

class PolynomialSplineContainer {
 public:

  using SplineType = PolynomialSplineQuintic;
  using SplineList = std::vector<SplineType>;

  PolynomialSplineContainer();
  virtual ~PolynomialSplineContainer();

  bool advance(double dt);
  bool addSpline(const SplineType& spline);
  bool addSpline(SplineType&& spline);
  bool reset();
  bool resetTime();

  double getContainerDuration() const;

  double getPosition() const;
  double getVelocity() const;
  double getAcceleration() const;

  double getPositionAtTime(double t) const;
  double getVelocityAtTime(double t) const;
  double getAccelerationAtTime(double t) const;

  double getEndPosition() const;
  double getEndVelocity() const;
  double getEndAcceleration() const;

  double getContainerTime() const;

  int getActiveSplineIndex() const;
  int getActiveSplineIndexAtTime(double t, double& timeOffset) const;
  bool isEmpty() const;

  virtual void setData(const std::vector<double>& knotPositions,
                       const std::vector<double>& knotValues,
                       double initialVelocity,
                       double initialAcceleration,
                       double finalVelocity,
                       double finalAcceleration);

  SplineType* getSpline(int splineIndex);

  void setContainerTime(double t);

  const SplineList& getSplines() const;

  static constexpr double undefinedValue = std::numeric_limits<double>::quiet_NaN();

 protected:
  int getCoeffIndex(int splineIdx, int aIdx) const;
  int getSplineColumnIndex(int splineIdx) const;

  SplineList splines_;
  double timeOffset_;
  double containerTime_;
  double containerDuration_;
  int activeSplineIdx_;
};

} /* namespace */
