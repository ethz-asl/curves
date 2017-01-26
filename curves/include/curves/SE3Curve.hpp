/*
 * ScalarCurveConfig.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Paul Furgale, Renaud Dube, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include "curves/SE3Config.hpp"
#include "curves/Curve.hpp"
#include <Eigen/Core>

namespace curves {

// Curves over SE3 inherit the interface from Curve and CurveBase and define specific
// methods to support physical interpretations of temporal derivatives.
//
// For the purposes of these uses, the curve is defined between Frame a and Frame b
// such that evaluate() returns \f$ \mathbf{T}_{a,b} \f$, the transformation that takes
// points from Frame b to Frame a.
//
class SE3Curve : public Curve<SE3Config> {
 public:
  SE3Curve();
  virtual ~SE3Curve();

  typedef Curve<SE3Config> Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;

  /// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularVelocityA(Time time) = 0;

  /// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularVelocityB(Time time) = 0;

  /// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearVelocityA(Time time) = 0;

  /// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearVelocityB(Time time) = 0;

  /// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistA(Time time) = 0;

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistB(Time time) = 0;

  /// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) = 0;

  /// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) = 0;

  /// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) = 0;

  /// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) = 0;

  /// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeA(unsigned derivativeOrder, Time time) = 0;

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeB(unsigned derivativeOrder, Time time) = 0;

  // Following functions added from SlerpSE3 curve .. todo clean

  /// \brief set the minimum sampling period
  virtual void setMinSamplingPeriod(Time time) = 0;

  /// \brief Set the sampling ratio.
  ///   eg. 4 will add a coefficient every 4 extend
  virtual void setSamplingRatio(const int ratio) = 0;

  virtual void clear() = 0;

  /// \brief Perform a rigid transformation on the left side of the curve
  virtual void transformCurve(const ValueType T) = 0;

  virtual void saveCurveTimesAndValues(const std::string& filename) const = 0;

  virtual void saveCurveAtTimes(const std::string& filename, std::vector<Time> times) const = 0;

  virtual void saveCorrectionCurveAtTimes(const std::string& filename, std::vector<Time> times) const = 0;

  virtual void getCurveTimes(std::vector<Time>* outTimes) const = 0;

  // Fake functions to comply with the current interfaces of trajectories_optimizer
  // todo : tidy up

  /// \brief Returns the number of coefficients in the correction curve
  virtual int correctionSize() const = 0;

  /// \brief Fold in the correction curve into the base curve and reinitialize
  ///        correction curve coefficients to identity transformations.
  virtual void foldInCorrections() = 0;

  /// \brief Add coefficients to the correction curve at given times.
  virtual void setCorrectionTimes(const std::vector<Time>& times) = 0;

  /// \brief Remove a correction coefficient at the specified time.
  virtual void removeCorrectionCoefficientAtTime(Time time) = 0;

  /// \brief Set the correction coefficient value at the specified time.
  virtual void setCorrectionCoefficientAtTime(Time time, ValueType value) = 0;

  /// \brief Reset the correction curve to identity values with knots at desired times
  virtual void resetCorrectionCurve(const std::vector<Time>& times) = 0;

  /// \brief Set the base curve to given values with knots at desired times
  /// Resets the curve beforehand.
  virtual void setBaseCurve(const std::vector<Time>& times, const std::vector<ValueType>& values) = 0;

  /// \brief Add / replace the given coefficients without resetting the curve.
  virtual void setBaseCurvePart(const std::vector<Time>& times, const std::vector<ValueType>& values) = 0;

  /// \brief Modifies values of the base coefficient in batch, starting at times[0] and assuming that
  /// a coefficient exists at all the specified times.
  virtual void modifyBaseCoefficientsValuesInBatch(const std::vector<Time>& times, const std::vector<ValueType>& values) = 0;

  virtual void getBaseCurveTimes(std::vector<Time>* outTimes) const = 0;

  virtual void getBaseCurveTimesInWindow(std::vector<Time>* outTimes, Time begTime, Time endTime) const = 0;

  virtual bool isEmpty() const = 0;

  virtual int size() const = 0;

  virtual int baseSize() const = 0;

  virtual void saveCorrectionCurveTimesAndValues(const std::string& filename) const = 0;

  /// \brief Get the dimension of this curve
  //virtual size_t dim() const;
  ///@}
 private:

};

} // namespace
