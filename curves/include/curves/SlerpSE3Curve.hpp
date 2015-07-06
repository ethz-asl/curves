/*
 * @file SlerpSE3Curve.hpp
 * @date Oct 10, 2014
 * @author Renaud Dube, Abel Gawel
 */

#ifndef CURVES_SLERP_SE3_CURVE_HPP
#define CURVES_SLERP_SE3_CURVE_HPP

#include "SE3Curve.hpp"
#include "LocalSupport2CoefficientManager.hpp"
#include "SE3CompositionCurve.hpp"
#include "gtsam/nonlinear/NonlinearFactorGraph.h"
#include "SamplingPolicy.hpp"

namespace curves {

/// Implements the Slerp (Spherical linear interpolation) curve class.
/// The Slerp interpolation function is defined as, with the respective Jacobians regarding  A and B:
/// \f[ T = A(A^{-1}B)^{\alpha} \f]
class SlerpSE3Curve : public SE3Curve {
  friend class SE3CompositionCurve<SlerpSE3Curve, SlerpSE3Curve>;
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

  virtual void clear();

  /// \brief Add factors to constrain the variables active at this time.
  void addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const;

  /// \brief Perform a rigid transformation on the left side of the curve
  void transformCurve(const ValueType T);

  virtual Time getTimeAtKey(gtsam::Key key) const;

 private:
  LocalSupport2CoefficientManager<Coefficient> manager_;
  SamplingPolicy slerpPolicy_;
};

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

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
      //todo: deal with extending curve with decreasing time
      // Get an iterator to the previous two coefficients
      SlerpSE3Curve::TimeToKeyCoefficientMap::iterator itPrev = (--curve->manager_.coefficientEnd());
      SlerpSE3Curve::TimeToKeyCoefficientMap::iterator itPrevPrev = (--(--curve->manager_.coefficientEnd()));

      Time tPrev = itPrev->first;
      Time tPrevPrev = itPrevPrev->first;
      if (tPrev - tPrevPrev >= minSamplingPeriod_) {
        // case 1 : the time delta between the two last knots is larger or equal to the minSamplingPeriod_
        // simply add a new coefficient and keep the previous one fixed
        curve->manager_.addCoefficientAtEnd(times[0], values[0], outKeys);
      } else if (times[0] - tPrevPrev > minSamplingPeriod_){
        //  add knot at tNew + move tPrev to tPrevPrev + minSamplingPeriod_ with value = interpolation
        curve->manager_.addCoefficientAtEnd(times[0], values[0], outKeys);
        SE3 newValue = curve->evaluate(tPrevPrev + minSamplingPeriod_);
        Time newTime = tPrevPrev + minSamplingPeriod_;
        curve->manager_.modifyCoefficient(itPrev, newTime, newValue);
      } else {
        // move knot at tNew with value = new value
        curve->manager_.modifyCoefficient(itPrev, times[0], values[0]);
      }
    }
  }
}

} // namespace curves

#endif /* CURVES_SLERP_SE3_CURVE_HPP */
