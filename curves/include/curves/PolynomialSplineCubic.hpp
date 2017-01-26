/*
 * PolynomialSplineCubic.hpp
 *
 *  Created on: Dec 20, 2014
 *      Author: dario, Peter Fankhauser
 */

#pragma once

#include "curves/PolynomialSplineBase.hpp"

namespace curves {

class PolynomialSplineCubic : public PolynomialSplineBase {
 public:
  PolynomialSplineCubic();
  virtual ~PolynomialSplineCubic();

  const std::vector<double>& getCoeffs() const;
  bool evalCoeffs(const SplineOpts& opts);
  void setCoeffsAndDuration(const std::vector<double>& coeffs, double duration);

  double getPositionAtTime(double dt) const;
  double getVelocityAtTime(double dt) const;
  double getAccelerationAtTime(double dt) const;

  void advanceTime(double dt);
  void resetTime();
  double getTime() const;

  double getSplineDuration() const;

 protected:
  double time_;
  double splineDuration_;
  bool didEvaluateCoeffs_;

  /*
  * s(t) = a3*t^3 + a2*t^2 + a1*t + a0
  * splineCoeff_ = [a0 a1 a2 a3]
  */
  std::vector<double> splineCoeff_;


};

} /* namespace */
