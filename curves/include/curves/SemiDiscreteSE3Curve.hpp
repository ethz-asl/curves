/*
 * @file SemiDiscreteSE3Curve.hpp
 * @date Feb 3, 2016
 * @author Renaud Dubé, Abel Gawel
 */

#ifndef CURVES_SEMI_DISCRETE_SE3_CURVE_HPP
#define CURVES_SEMI_DISCRETE_SE3_CURVE_HPP

#include "SE3Curve.hpp"
#include "LocalSupport2CoefficientManager.hpp"
#include "SE3CompositionCurve.hpp"
#include "gtsam/nonlinear/NonlinearFactorGraph.h"
#include "SamplingPolicy.hpp"
#include "CubicHermiteSE3Curve.hpp"

namespace curves {

/// Implements a discrete SE3 curve class.
class SemiDiscreteSE3Curve : public SE3Curve {
  friend class SE3CompositionCurve<SemiDiscreteSE3Curve, SemiDiscreteSE3Curve>;
  friend class SamplingPolicy;
 public:
  typedef ValueType Coefficient;
  typedef LocalSupport2CoefficientManager<Coefficient>::TimeToKeyCoefficientMap TimeToKeyCoefficientMap;
  typedef LocalSupport2CoefficientManager<Coefficient>::CoefficientIter CoefficientIter;

  SemiDiscreteSE3Curve();
  virtual ~SemiDiscreteSE3Curve();

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

  /// \brief Get an evaluator at this time
  virtual gtsam::Expression<ValueType> getValueExpression(const Time& time) const;

  virtual gtsam::Expression<DerivativeType> getDerivativeExpression(const Time& time, unsigned derivativeOrder) const;

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

  /// Initialize a GTSAM values structure with the desired keys
  virtual void initializeGTSAMValues(gtsam::KeySet keys, gtsam::Values* values) const;

  /// Initialize a GTSAM values structure for all keys
  virtual void initializeGTSAMValues(gtsam::Values* values) const;

  // updates the relevant curve coefficients from the GTSAM values structure
  virtual void updateFromGTSAMValues(const gtsam::Values& values);

  // set minimum sampling period
  void setMinSamplingPeriod(Time time);

  /// \brief Set the sampling ratio.
  ///   eg. 4 will add a coefficient every 4 extend
  void setSamplingRatio(const int ratio);

  virtual void clear();

  /// \brief Add factors to constrain the variables active at this time.
  void addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const;

  /// \brief Perform a rigid transformation on the left side of the curve
  void transformCurve(const ValueType T);

  virtual Time getTimeAtKey(gtsam::Key key) const;

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
  SamplingPolicy discretePolicy_;
};

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;


// extend policy for slerp curves
template<>
inline void SamplingPolicy::extend<SemiDiscreteSE3Curve, SE3>(const std::vector<Time>& times,
                                                          const std::vector<SE3>& values,
                                                          SemiDiscreteSE3Curve* curve,
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
          SemiDiscreteSE3Curve::TimeToKeyCoefficientMap::iterator itPrev = (--curve->manager_.coefficientEnd());
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

#endif /* CURVES_DISCRETE_SE3_CURVE_HPP */
