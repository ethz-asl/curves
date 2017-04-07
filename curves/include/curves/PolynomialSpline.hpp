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
#include "curves/polynomial_splines_traits.hpp"

namespace curves {

template <int splineOrder_>
class PolynomialSpline {
 public:

  static constexpr unsigned int splineOrder = splineOrder_;
  static constexpr unsigned int coefficientCount = splineOrder + 1;

  using SplineImplementation = spline_traits::spline_rep<double, splineOrder>;
  using SplineCoefficients = typename SplineImplementation::SplineCoefficients;
  using EigenTimeVectorType = Eigen::Matrix<double, 1, coefficientCount>;
  using EigenCoefficientVectorType = Eigen::Matrix<double, coefficientCount, 1>;

  PolynomialSpline() :
    currentTime_(0.0),
    duration_(0.0),
    didEvaluateCoeffs_(false),
    coefficients_()
  {

  }

  virtual ~PolynomialSpline() {

  }

  const SplineCoefficients& getCoefficients() const {
    return coefficients_;
  }

  bool computeCoefficients(const SplineOptions& options) {
    SplineImplementation::compute(options, coefficients_);
    duration_ = options.tf_;
    return true;
  }

  void setCoefficientsAndDuration(const SplineCoefficients& coefficients, double duration) {
    coefficients_ = coefficients;
    duration_ = duration;
  }

  void setCoefficientsAndDuration(const EigenCoefficientVectorType& coefficients, double duration) {
    for (unsigned int k=0; k<coefficientCount; k++) {
      coefficients_[k] = coefficients(k);
    }
    duration_ = duration;
  }

  constexpr double getPositionAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(),
                              SplineImplementation::tau(std::max(0.0, std::min(tk, duration_))).begin(), 0.0);
  }

  constexpr double getVelocityAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(),
                              SplineImplementation::dtau(std::max(0.0, std::min(tk, duration_))).begin(), 0.0);
  }

  constexpr double getAccelerationAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(),
                              SplineImplementation::ddtau(std::max(0.0, std::min(tk, duration_))).begin(), 0.0);
  }

  static inline void getTimeVector(Eigen::Ref<EigenTimeVectorType> timeVec, double tk) {
    timeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::tau(tk).data());
  }

  static inline void getdTimeVector(Eigen::Ref<EigenTimeVectorType> dtimeVec, double tk) {
    dtimeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::dtau(tk).data());
  }

  static inline void getddTimeVector(Eigen::Ref<EigenTimeVectorType> ddtimeVec, double tk) {
    ddtimeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::ddtau(tk).data());
  }

  static inline void getTimeVectorAtZero(Eigen::Ref<EigenTimeVectorType> timeVec) {
    timeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::tauZero.data());
  }

  static inline void getdTimeVectorAtZero(Eigen::Ref<EigenTimeVectorType> dtimeVec) {
    dtimeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::dtauZero.data());
  }

  static inline void getddTimeVectorAtZero(Eigen::Ref<EigenTimeVectorType> ddtimeVec) {
    ddtimeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::ddtauZero.data());
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
   * s(t) = an*t^n + ... + a1*t + a0
   * splineCoeff_ = [an ... a1 a0]
   */
  SplineCoefficients coefficients_;
};

} /* namespace */
