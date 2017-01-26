/*
 * SlerpSE3Curve.hpp
 *
 *  Created on: Oct 10, 2014
 *      Author: Renaud Dube, Abel Gawel, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include "SE3Curve.hpp"
#include "LocalSupport2CoefficientManager.hpp"
#include "kindr/Core"
#include "SE3CompositionCurve.hpp"
#include "SamplingPolicy.hpp"
#include "CubicHermiteSE3Curve.hpp"

namespace curves {

/// Implements the Slerp (Spherical linear interpolation) curve class.
/// The Slerp interpolation function is defined as, with the respective Jacobians regarding  A and B:
/// \f[ T = A(A^{-1}B)^{\alpha} \f]
class SlerpSE3Curve : public SE3Curve {
  friend class SE3CompositionCurve<SlerpSE3Curve, SlerpSE3Curve>;
  friend class SE3CompositionCurve<SlerpSE3Curve, CubicHermiteSE3Curve>;
  friend class SamplingPolicy;
 public:
  typedef SE3Curve::ValueType ValueType;
  typedef SE3Curve::DerivativeType DerivativeType;
  typedef ValueType Coefficient;
  typedef LocalSupport2CoefficientManager<Coefficient>::TimeToKeyCoefficientMap TimeToKeyCoefficientMap;
  typedef LocalSupport2CoefficientManager<Coefficient>::CoefficientIter CoefficientIter;

  SlerpSE3Curve();
  virtual ~SlerpSE3Curve();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// The first valid time for the curve.
  virtual Time getMinTime() const;

  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const;

  bool isEmpty() const;

  // return number of coefficients curve is composed of
  int size() const;

  /// Extend the curve so that it can be evaluated at these times.
  /// Try to make the curve fit to the values.
  /// Underneath the curve should have some default policy for fitting.
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys = NULL);

  /// \brief Fit a new curve to these data points.
  ///
  /// The existing curve will be cleared.
  /// Underneath the curve should have some default policy for fitting.
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL);

  /// \brief Set some coefficients of the curve
  /// The existing curve will NOT be cleared.
  void setCurve(const std::vector<Time>& times,
                const std::vector<ValueType>& values);

  /// Evaluate the ambient space of the curve.
  virtual ValueType evaluate(Time time) const;

  /// Evaluate the curve derivatives.
  /// linear 1st derivative has following behaviour:
  /// - time is out of bound --> error
  /// - time is between 2 coefficients --> take slope between the 2 coefficients
  /// - time is on coefficient (not last coefficient) --> take slope between coefficient and next coefficients
  /// - time is on last coefficient --> take slope between last-1 and last coefficient
  /// derivatives of order >1 equal 0
  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) const;

  virtual void setTimeRange(Time minTime, Time maxTime);

  /// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularVelocityA(Time time);

  /// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularVelocityB(Time time);

  /// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearVelocityA(Time time);

  /// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearVelocityB(Time time);

  /// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistA(Time time);

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistB(Time time);

  /// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularDerivativeB(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearDerivativeB(unsigned derivativeOrder, Time time);

  /// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeB(unsigned derivativeOrder, Time time);

  // set minimum sampling period
  void setMinSamplingPeriod(Time time);

  /// \brief Set the sampling ratio.
  ///   eg. 4 will add a coefficient every 4 extend
  void setSamplingRatio(const int ratio);

  virtual void clear();

  /// \brief Perform a rigid transformation on the left side of the curve
  void transformCurve(const ValueType T);

  void saveCurveTimesAndValues(const std::string& filename) const;

  void saveCurveAtTimes(const std::string& filename, std::vector<Time> times) const;

  void saveCorrectionCurveAtTimes(const std::string& filename, std::vector<Time> times) const {};

  void getCurveTimes(std::vector<Time>* outTimes) const;

  // Fake functions to comply with the current interfaces of trajectories_optimizer
  // todo : tidy up

  /// \brief Returns the number of coefficients in the correction curve
  int correctionSize() const {return 0;};

  /// \brief Fold in the correction curve into the base curve and reinitialize
  ///        correction curve coefficients to identity transformations.
  void foldInCorrections() {};

  /// \brief Add coefficients to the correction curve at given times.
  void setCorrectionTimes(const std::vector<Time>& times) {};

  /// \brief Remove a correction coefficient at the specified time.
  void removeCorrectionCoefficientAtTime(Time time) {};

  /// \brief Set the correction coefficient value at the specified time.
  void setCorrectionCoefficientAtTime(Time time, ValueType value) {};

  /// \brief Reset the correction curve to identity values with knots at desired times
  void resetCorrectionCurve(const std::vector<Time>& times) {};

  /// \brief Set the base curve to given values with knots at desired times
  /// Resets the curve beforehand.
  void setBaseCurve(const std::vector<Time>& times, const std::vector<ValueType>& values) {};

  /// \brief Add / replace the given coefficients without resetting the curve.
  void setBaseCurvePart(const std::vector<Time>& times, const std::vector<ValueType>& values) {};

  /// \brief Modifies values of the base coefficient in batch, starting at times[0] and assuming that
  /// a coefficient exists at all the specified times.
  void modifyBaseCoefficientsValuesInBatch(const std::vector<Time>& times, const std::vector<ValueType>& values) {};

  void getBaseCurveTimes(std::vector<Time>* outTimes) const {};

  void getBaseCurveTimesInWindow(std::vector<Time>* outTimes, Time begTime, Time endTime) const {};

  // return number of coefficients curve is composed of
  int baseSize() const {return size();};

  void saveCorrectionCurveTimesAndValues(const std::string& filename) const {};

 private:
  LocalSupport2CoefficientManager<Coefficient> manager_;
  SamplingPolicy slerpPolicy_;
};

typedef kindr::HomogeneousTransformationPosition3RotationQuaternionD SE3;
typedef SE3::Rotation SO3;
typedef kindr::AngleAxisPD AngleAxis;

SE3 transformationPower(SE3  T, double alpha);

SE3 composeTransformations(SE3 A, SE3 B);

SE3 inverseTransformation(SE3 T);

// extend policy for slerp curves
template<>
inline void SamplingPolicy::extend<SlerpSE3Curve, SE3>(const std::vector<Time>& times,
                                                       const std::vector<SE3>& values,
                                                       SlerpSE3Curve* curve,
                                                       std::vector<Key>* outKeys) {
  //todo: deal with minSamplingPeriod_ when extending with multiple times
  if (times.size() != 1) {
    curve->manager_.insertCoefficients(times, values, outKeys);
  } else {
    //If the curve is empty or of size 1, simply add the new coefficient
    if (curve->isEmpty() || curve->size() == 1) {
      curve->manager_.insertCoefficients(times, values, outKeys);
    } else {
      if (minimumMeasurements_ == 1) {
        curve->manager_.addCoefficientAtEnd(times[0], values[0], outKeys);
      } else {
        ++measurementsSinceLastExtend_;

        if (measurementsSinceLastExtend_ == 1) {
          curve->manager_.addCoefficientAtEnd(times[0], values[0], outKeys);
        } else {
          SlerpSE3Curve::TimeToKeyCoefficientMap::iterator itPrev = (--curve->manager_.coefficientEnd());
          curve->manager_.modifyCoefficient(itPrev, times[0], values[0]);
        }
        if (measurementsSinceLastExtend_ == minimumMeasurements_) {
          measurementsSinceLastExtend_ = 0;
        }
      }
    }
  }
}

} // namespace curves

