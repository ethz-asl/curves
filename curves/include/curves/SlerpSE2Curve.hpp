/*
 * @file SlerpSE2Curve.hpp
 * @date Nov 24, 2015
 * @author Renaud Dubé, Abel Gawel
 */

#ifndef CURVES_SLERP_SE2_CURVE_HPP
#define CURVES_SLERP_SE2_CURVE_HPP

#include "SE2Curve.hpp"
#include "LocalSupport2CoefficientManager.hpp"
//#include "SE2CompositionCurve.hpp"
#include "gtsam/nonlinear/NonlinearFactorGraph.h"
#include "SamplingPolicy.hpp"

namespace curves {

/// Implements the Slerp (Spherical linear interpolation) curve class.
/// The Slerp interpolation function is defined as, with the respective Jacobians regarding  A and B:
/// \f[ T = A(A^{-1}B)^{\alpha} \f]
class SlerpSE2Curve : public SE2Curve {
  //friend class SE2CompositionCurve<SlerpSE2Curve, SlerpSE2Curve>;
  //friend class SE2CompositionCurve<SlerpSE2Curve, CubicHermiteSE2Curve>;
  friend class SamplingPolicy;
 public:
  typedef SE2Curve::ValueType ValueType;
  typedef SE2Curve::DerivativeType DerivativeType;
  typedef ValueType Coefficient;
  typedef LocalSupport2CoefficientManager<Coefficient>::TimeToKeyCoefficientMap TimeToKeyCoefficientMap;
  typedef LocalSupport2CoefficientManager<Coefficient>::CoefficientIter CoefficientIter;

  SlerpSE2Curve();
  virtual ~SlerpSE2Curve();

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

 private:
  LocalSupport2CoefficientManager<Coefficient> manager_;
  SamplingPolicy slerpPolicy_;
};

typedef gtsam::Pose2 SE2;
typedef gtsam::Rot2 SO2;

SE2 transformationPower(SE2  T, double alpha);

SE2 composeTransformations(SE2 A, SE2 B);

SE2 inverseTransformation(SE2 T);

// extend policy for slerp curves
template<>
inline void SamplingPolicy::extend<SlerpSE2Curve, SE2>(const std::vector<Time>& times,
                                                const std::vector<SE2>& values,
                                                SlerpSE2Curve* curve,
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
          SlerpSE2Curve::TimeToKeyCoefficientMap::iterator itPrev = (--curve->manager_.coefficientEnd());
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

#endif /* CURVES_SLERP_SE2_CURVE_HPP */
