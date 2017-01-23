/*
 * PolynomialSplineBase.hpp
 *
 *  Created on: Dec 1, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#pragma once

#include <vector>

namespace curves {

class PolynomialSplineBase {
 public:

  struct SplineOpts {
    double tf;
    double pos0 = 0.0;
    double posT = 1.0;
    double vel0 = 0.0;
    double velT = 0.0;
    double acc0 = 0.0;
    double accT = 0.0;
  };

  PolynomialSplineBase();
  virtual ~PolynomialSplineBase();

  virtual const std::vector<double>& getCoeffs() const = 0;
  virtual bool evalCoeffs(const SplineOpts& opts) = 0;
  virtual void setCoeffsAndDuration(const std::vector<double>& coeffs_, double duration) = 0;

  virtual double getPositionAtTime(double dt) const = 0;
  virtual double getVelocityAtTime(double dt) const = 0;
  virtual double getAccelerationAtTime(double dt) const = 0;

  virtual double getSplineDuration() const = 0;

};

} /* namespace */

