/*
 * PolynomialSpline.hpp
 *
 *  Created on: Mar 7, 2017
 *      Author: Dario Bellicoso
 */

#pragma once

// eigen
#include <Eigen/Core>

// curves
#include "curves/polynomial_splines.hpp"
#include "curves/polynomial_splines_traits.hpp"

namespace curves {

template <unsigned int splineOrder_>
class PolynomialSpline {
 public:

  constexpr unsigned int splineOrder = splineOrder_;
  constexpr unsigned int coefficientCount = splineOrder + 1;

  using SplineCoefficients = std::vector<double, coefficientCount>;

  using TimeVector = spline_traits::time_vector<double, this->splineOrder>;

  PolynomialSpline();
  virtual ~PolynomialSpline();

  const SplineCoefficients& getCoeffs() const {
    return coefficients_;
  }

  bool computeCoefficients(const SplineOptions& options) {
    return true;
  }

  void setCoefficientsAndDuration(const SplineCoefficients& coefficients, double duration) {
    coefficients_ = coefficients;
    duration_ = duration;
  }

  void setCoefficientsAndDuration(const Eigen::Matrix<double, coefficientCount, 1>& coefficients, double duration) {
    for (unsigned int k=0; k<coefficientCount; k++) {
      coefficients_[k] = coefficients(k);
    }
    duration_ = duration;
  }

  double getPositionAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(), TimeVector::tau(tk).begin(), 0.0);
  }

  double getVelocityAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(), TimeVector::dtau(tk).begin(), 0.0);
  }

  double getAccelerationAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(), TimeVector::ddtau(tk).begin(), 0.0);
  }

  void advanceTime(double dt) {
    currentTime_ += dt;
  }

  void resetTime() {
    currentTime_ = 0.0;
  }

  double getTime() const {
    return currentTime_;
  }

  double getSplineDuration() const {
    return duration_;
  }

 protected:
  //! A helper counter which used to get the state of the spline.
  double currentTime_;

  //! The duration of the spline in seconds.
  double duration_;

  //! True if the coefficents were computed at least once.
  bool didEvaluateCoeffs_;

  /*
   * s(t) = a5*t^5 + a4*t^5 + a3*t^3 + a2*t^2 + a1*t + a0
   * splineCoeff_ = [a0 a1 a2 a3 a4 a5]
   */
  SplineCoefficients coefficients_;
};

using PolynomialSplineThird = PolynomialSpline<3>;
using PolynomialSplineFifth = PolynomialSpline<5>;

} /* namespace */
