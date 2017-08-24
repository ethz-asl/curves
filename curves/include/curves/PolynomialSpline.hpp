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

/*
 * This class is the implementation of a scalar polynomial spline s(t) function of a scalar t.
 * The spline is define as
 *    s(t) = an*t^n + ... + a1*t + a0 = sum(ai*t^i)
 *
 *  The spline coefficients are stored in a standard container as
 *    alpha = [an ... a1 a0]
 */
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
    duration_(0.0),
    didEvaluateCoeffs_(false),
    coefficients_()
  {

  }

  template<typename SplineCoeff_>
  PolynomialSpline(SplineCoeff_&& coefficients, double duration) :
    duration_(duration),
    didEvaluateCoeffs_(true),
    coefficients_(std::forward<SplineCoeff_>(coefficients))
  {

  }

  explicit PolynomialSpline(const SplineOptions& options) : duration_(options.tf_) {
    computeCoefficients(options);
  }

  explicit PolynomialSpline(SplineOptions&& options) : duration_(options.tf_) {
    computeCoefficients(std::move(options));
  }

  virtual ~PolynomialSpline() { }

  PolynomialSpline(PolynomialSpline &&) = default;
  PolynomialSpline& operator=(PolynomialSpline &&) = default;

  PolynomialSpline(const PolynomialSpline&) = default;
  PolynomialSpline& operator=(const PolynomialSpline&) = default;

  //! Get the coefficients of the spline.
  const SplineCoefficients& getCoefficients() const {
    return coefficients_;
  }

  //! Compute the coefficients of the spline.
  template<typename SplineOptionsType_>
  bool computeCoefficients(SplineOptionsType_&& options) {
    duration_ = options.tf_;
    return SplineImplementation::compute(std::forward<SplineOptionsType_>(options), coefficients_);
  }

  //! Set the coefficients and the duration of the spline.
  void setCoefficientsAndDuration(const SplineCoefficients& coefficients, double duration) {
    coefficients_ = coefficients;
    duration_ = duration;
  }

  //! Get the spline evaluated at time tk.
  constexpr double getPositionAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(),
                              SplineImplementation::tau(std::max(0.0, std::min(tk, duration_))).begin(), 0.0);
  }

  //! Get the first derivative of the spline evaluated at time tk.
  constexpr double getVelocityAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(),
                              SplineImplementation::dtau(std::max(0.0, std::min(tk, duration_))).begin(), 0.0);
  }

  //! Get the second derivative of the spline evaluated at time tk.
  constexpr double getAccelerationAtTime(double tk) const {
    return std::inner_product(coefficients_.begin(), coefficients_.end(),
                              SplineImplementation::ddtau(std::max(0.0, std::min(tk, duration_))).begin(), 0.0);
  }

  //! Get the time vector tau evaluated at time tk.
  static inline void getTimeVector(Eigen::Ref<EigenTimeVectorType> timeVec, double tk) {
    timeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::tau(tk).data());
  }

  //! Get the first derivative of the time vector tau evaluated at time tk.
  static inline void getdTimeVector(Eigen::Ref<EigenTimeVectorType> dtimeVec, double tk) {
    dtimeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::dtau(tk).data());
  }

  //! Get the second derivative of the time vector tau evaluated at time tk.
  static inline void getddTimeVector(Eigen::Ref<EigenTimeVectorType> ddtimeVec, double tk) {
    ddtimeVec = Eigen::Map<EigenTimeVectorType>(SplineImplementation::ddtau(tk).data());
  }

  //! Get the time vector tau evaluated at zero.
  static inline void getTimeVectorAtZero(Eigen::Ref<EigenTimeVectorType> timeVec) {
    timeVec = Eigen::Map<const EigenTimeVectorType>((SplineImplementation::tauZero).data());
  }

  //! Get the first derivative of the time vector tau evaluated at zero.
  static inline void getdTimeVectorAtZero(Eigen::Ref<EigenTimeVectorType> dtimeVec) {
    dtimeVec = Eigen::Map<const EigenTimeVectorType>((SplineImplementation::dtauZero).data());
  }

  //! Get the second derivative of the time vector tau evaluated at zero.
  static inline void getddTimeVectorAtZero(Eigen::Ref<EigenTimeVectorType> ddtimeVec) {
    ddtimeVec = Eigen::Map<const EigenTimeVectorType>((SplineImplementation::ddtauZero).data());
  }

  //! Get the duration of the spline in seconds.
  double getSplineDuration() const {
    return duration_;
  }

 protected:
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
