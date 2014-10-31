#ifndef CURVES_TYPED_CURVE_HPP
#define CURVES_TYPED_CURVE_HPP

#include "CurveBase.hpp"
#include "Evaluator.hpp"

#include "gtsam_unstable/nonlinear/Expression.h"

namespace curves {

template<typename CurveConfig>
class Curve : public CurveBase
{
 public:
  
  /// The value type of the curve.
  typedef typename CurveConfig::ValueType ValueType;

  /// The curve's derivative type.
  typedef typename CurveConfig::DerivativeType DerivativeType;

  typedef Evaluator<CurveConfig> EvaluatorType;
  
  typedef typename EvaluatorType::Ptr EvaluatorTypePtr;

  Curve() { }
  virtual ~Curve() { }

  /// \name Methods to evaluate the curve
  ///@{

  /// Evaluate the ambient space of the curve.
  virtual ValueType evaluate(Time time) const = 0;
  
  /// Evaluate the curve derivatives.
  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) const = 0;

//  /// \brief Get an evaluator at this time.
//  virtual EvaluatorTypePtr getEvaluator(const Time& time) const = 0;

  /// \brief Get a gtsam::Expression which evaluates the curve at this time.
  virtual gtsam::Expression<ValueType> getEvalExpression(const Time& time) const = 0;

  ///@}

  /// \name Methods to fit the curve based on data.
  ///@{

  /// Extend the curve so that it can be evaluated at these times.
  /// Try to make the curve fit to the values.
  /// Underneath the curve should have some default policy for fitting.
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values) = 0;

  /// \brief Fit a new curve to these data points.
  ///
  /// The existing curve will be cleared.
  /// Underneath the curve should have some default policy for fitting.
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values) = 0;

  ///@}

};

} // namespace curves


#endif /* CURVES_TYPED_CURVE_HPP */
