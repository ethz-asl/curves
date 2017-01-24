/*
 * PolynomialSplineCubic.cpp
 *
 *  Created on: Dec 20, 2014
 *      Author: dario, Peter Fankhauser
 */

#include "curves/PolynomialSplineCubic.hpp"
#include <Eigen/Dense>

namespace curves {

PolynomialSplineCubic::PolynomialSplineCubic():
      time_(0.0),
      splineDuration_(0.0),
      didEvaluateCoeffs_(false),
      splineCoeff_(3)
{
}


PolynomialSplineCubic::~PolynomialSplineCubic() {
}


const std::vector<double>& PolynomialSplineCubic::getCoeffs() const {
  return splineCoeff_;
}


bool PolynomialSplineCubic::evalCoeffs(const SplineOpts& opts) {
  didEvaluateCoeffs_ = false;

  Eigen::Matrix<double,4,1> b;
  b << opts.pos0, 0.0, opts.posT, 0.0;

  Eigen::Matrix<double,6,6> A;
  A << 1.0,     0.0,        0.0,              0.0,                  0.0,                    0.0,
       0.0,     1.0,        0.0,              0.0,                  0.0,                    0.0,
       0.0,     0.0,        2.0,              0.0,                  0.0,                    0.0,
       1.0,     opts.tf,    pow(opts.tf,2),   pow(opts.tf,3),       pow(opts.tf,4),         pow(opts.tf,5),
       0.0,     1.0,        2.0*opts.tf,      3.0*pow(opts.tf,2),   4.0*pow(opts.tf,3),     5.0*pow(opts.tf,4),
       0.0,     0.0,        2.0,              6.0*opts.tf,          12.0*pow(opts.tf,2),    20.0*pow(opts.tf,3);

  Eigen::VectorXd coeffs;
  coeffs = A.colPivHouseholderQr().solve(b);

  // copy result in std::vector array
  for (int k=0; k<4; k++) {
    splineCoeff_[k] = coeffs(k);
  }


  // save spline options
  splineDuration_ = opts.tf;

  didEvaluateCoeffs_ = true;

  return didEvaluateCoeffs_;
}


void PolynomialSplineCubic::setCoeffsAndDuration(const std::vector<double>& coeffs, double duration) {
  splineCoeff_ = coeffs;
  splineDuration_ = duration;
}

double PolynomialSplineCubic::getPositionAtTime(double dt) const {
  return splineCoeff_[3]*pow(dt,3) + splineCoeff_[2]*pow(dt,2) + splineCoeff_[1]*dt + splineCoeff_[0];
}

double PolynomialSplineCubic::getVelocityAtTime(double dt) const {
  return 3.0*splineCoeff_[3]*pow(dt,2) + 2.0*splineCoeff_[2]*dt + splineCoeff_[1];
}

double PolynomialSplineCubic::getAccelerationAtTime(double dt) const {
  return 6.0*splineCoeff_[3]*(dt) + 2.0*splineCoeff_[2];
}


void PolynomialSplineCubic::advanceTime(double dt) {
  time_ += dt;
}


void PolynomialSplineCubic::resetTime() {
  time_ = 0.0;
}


double PolynomialSplineCubic::getTime() const {
  return time_;
}

double PolynomialSplineCubic::getSplineDuration() const {
  return splineDuration_;
}

} /* namespace */
