/*
 * PolynomialSplineQuintic.cpp
 *
 *  Created on: Dec 1, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#include "curves/PolynomialSplineQuintic.hpp"
#include <Eigen/Dense>

// boost
#include <boost/math/special_functions/pow.hpp>


static const int NUM_SPLINE_COEFFS = 6;

namespace curves {

PolynomialSplineQuintic::PolynomialSplineQuintic():
    time_(0.0),
    splineDuration_(0.0),
    didEvaluateCoeffs_(false),
    splineCoeff_(NUM_SPLINE_COEFFS, 0.0)
{
}

PolynomialSplineQuintic::~PolynomialSplineQuintic() {
}

void PolynomialSplineQuintic::advanceTime(double dt) {
  time_ += dt;
}

void PolynomialSplineQuintic::resetTime() {
  time_ = 0.0;
}

double PolynomialSplineQuintic::getTime() const {
  return time_;
}

const PolynomialSplineQuintic::SplineCoefficients& PolynomialSplineQuintic::getCoeffs() const {
  return splineCoeff_;
}

bool PolynomialSplineQuintic::evalCoeffs(const SplineOpts& opts) {
  using namespace boost::math;

  didEvaluateCoeffs_ = false;

  Eigen::Matrix<double, NUM_SPLINE_COEFFS, 1> b;
  b << opts.pos0, opts.vel0, opts.acc0, opts.posT, opts.velT, opts.accT;

  Eigen::Matrix<double,6,6> A;
  A << 1.0,     0.0,        0.0,              0.0,                  0.0,                    0.0,
       0.0,     1.0,        0.0,              0.0,                  0.0,                    0.0,
       0.0,     0.0,        2.0,              0.0,                  0.0,                    0.0,
       1.0,     opts.tf,    pow<2>(opts.tf),  pow<3>(opts.tf),      pow<4>(opts.tf),        pow<5>(opts.tf),
       0.0,     1.0,        2.0*opts.tf,      3.0*pow<2>(opts.tf),  4.0*pow<3>(opts.tf),    5.0*pow<4>(opts.tf),
       0.0,     0.0,        2.0,              6.0*opts.tf,          12.0*pow<2>(opts.tf),   20.0*pow<3>(opts.tf);

  Eigen::VectorXd coeffs;
  coeffs = A.colPivHouseholderQr().solve(b);

  // copy result in std::vector array
  for (int k=0; k<NUM_SPLINE_COEFFS; k++) {
    splineCoeff_[k] = coeffs(k);
  }

  // save spline options
  splineDuration_ = opts.tf;

  didEvaluateCoeffs_ = true;

  return didEvaluateCoeffs_;
}

void PolynomialSplineQuintic::setCoeffsAndDuration(const PolynomialSplineQuintic::SplineCoefficients& coeffs, double duration) {
  splineCoeff_ = coeffs;
  splineDuration_ = duration;
}

void PolynomialSplineQuintic::setCoeffsAndDuration(const Eigen::Matrix<double, 6, 1>& coeffs, double duration) {
  for (unsigned int k=0; k<6; k++) {
    splineCoeff_[k] = coeffs(k);
  }
  splineDuration_ = duration;
}


double PolynomialSplineQuintic::getPositionAtTime(double tk) const {

  tk = std::max(0.0, std::min(tk, splineDuration_));

  // Direct multiplication is more efficient then pow(x, y).
//  return   splineCoeff_[5]*boost::math::pow<5>(tk)
//         + splineCoeff_[4]*boost::math::pow<4>(tk)
//         + splineCoeff_[3]*boost::math::pow<3>(tk)
//         + splineCoeff_[2]*tk*tk
//         + splineCoeff_[1]*tk
//         + splineCoeff_[0];
  return   splineCoeff_[5]*tk*tk*tk*tk*tk
         + splineCoeff_[4]*tk*tk*tk*tk
         + splineCoeff_[3]*tk*tk*tk
         + splineCoeff_[2]*tk*tk
         + splineCoeff_[1]*tk
         + splineCoeff_[0];
}

double PolynomialSplineQuintic::getVelocityAtTime(double tk) const {

  tk = std::max(0.0, std::min(tk, splineDuration_));

//  return   5.0*splineCoeff_[5]*boost::math::pow<4>(tk)
//         + 4.0*splineCoeff_[4]*boost::math::pow<3>(tk)
//         + 3.0*splineCoeff_[3]*tk*tk
//         + 2.0*splineCoeff_[2]*tk
//         + splineCoeff_[1];

  return   5.0*splineCoeff_[5]*tk*tk*tk*tk
         + 4.0*splineCoeff_[4]*tk*tk*tk
         + 3.0*splineCoeff_[3]*tk*tk
         + 2.0*splineCoeff_[2]*tk
         + splineCoeff_[1];
}

double PolynomialSplineQuintic::getAccelerationAtTime(double tk) const {

  tk = std::max ( 0.0, std::min(tk, splineDuration_) );

//  return   20.0*splineCoeff_[5]*boost::math::pow<3>(tk)
//         + 12.0*splineCoeff_[4]*tk*tk
//         +  6.0*splineCoeff_[3]*tk
//         +  2.0*splineCoeff_[2];

  return   20.0*splineCoeff_[5]*tk*tk*tk
         + 12.0*splineCoeff_[4]*tk*tk
         +  6.0*splineCoeff_[3]*tk
         +  2.0*splineCoeff_[2];
}

double PolynomialSplineQuintic::getSplineDuration() const {
  return splineDuration_;
}


} /* namespace */
