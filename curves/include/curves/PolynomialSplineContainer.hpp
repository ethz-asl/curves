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
  bool addSpline(const PolynomialSplineQuintic& spline);
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

  PolynomialSplineBase* getSpline(int splineIndex);

  void setContainerTime(double t);

  const std::vector<PolynomialSplineQuintic>& getSplines() const;

  static constexpr double undefinedValue = std::numeric_limits<double>::quiet_NaN();

 protected:
  void getTimeVector(Eigen::Matrix<double, 1, 6>& timeVec, double t_k) const;
  void getdTimeVector(Eigen::Matrix<double, 1, 6>& timeVec, double t_k) const;
  void getddTimeVector(Eigen::Matrix<double, 1, 6>& timeVec, double t_k) const;

  int getCoeffIndex(int splineIdx, int aIdx) const;
  int getSplineColumnIndex(int splineIdx) const;

  std::vector<PolynomialSplineQuintic> splines_;
  double timeOffset_;
  double containerTime_;
  double containerDuration_;
  int activeSplineIdx_;
};

} /* namespace */
