/*
 * PolynomialSplineQuintic.hpp
 *
 *  Created on: Dec 1, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#pragma once

#include "curves/PolynomialSplineBase.hpp"
#include <Eigen/Core>

namespace curves {

class PolynomialSplineQuintic : public PolynomialSplineBase {
 public:

  using SplineCoefficients = std::vector<double>;

  PolynomialSplineQuintic();
  virtual ~PolynomialSplineQuintic();

  const SplineCoefficients& getCoeffs() const;
  bool evalCoeffs(const SplineOpts& opts);
  void setCoeffsAndDuration(const SplineCoefficients& coeffs, double duration);
  void setCoeffsAndDuration(const Eigen::Matrix<double, 6, 1>& coeffs, double duration);

  double getPositionAtTime(double tk) const;
  double getVelocityAtTime(double tk) const;
  double getAccelerationAtTime(double tk) const;

  void advanceTime(double dt);
  void resetTime();
  double getTime() const;

  double getSplineDuration() const;

 protected:
  double time_;
  double splineDuration_;
  bool didEvaluateCoeffs_;

  /*
   * s(t) = a5*t^5 + a4*t^5 + a3*t^3 + a2*t^2 + a1*t + a0
   * splineCoeff_ = [a0 a1 a2 a3 a4 a5]
   */
  SplineCoefficients splineCoeff_;
};

} /* namespace */
