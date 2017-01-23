/*
 * PolynomialSplineContainer.hpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#pragma once

#include "curves/PolynomialSplineQuintic.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits>

namespace curves {

class PolynomialSplineContainer {
 public:
  PolynomialSplineContainer();
  virtual ~PolynomialSplineContainer();

  bool advance(double dt);
  bool addSpline(PolynomialSplineQuintic& spline);
  bool reset();
  bool resetTime();

  double getContainerDuration() const;

  double getPosition();
  double getVelocity();
  double getAcceleration();

  double getPositionAtTime(double t) const;
  double getVelocityAtTime(double t) const;
  double getAccelerationAtTime(double t) const;

  double getEndPosition();
  double getEndVelocity();
  double getEndAcceleration();

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

  PolynomialSplineBase* getSpline(int splineIndex);

  void setContainerTime(double t);

  static constexpr double undefinedValue = std::numeric_limits<double>::quiet_NaN();

 protected:
  std::vector<PolynomialSplineQuintic> splines_;
  double timeOffset_;
  double containerTime_;
  double containerDuration_;
  int activeSplineIdx_;
};

} /* namespace */
