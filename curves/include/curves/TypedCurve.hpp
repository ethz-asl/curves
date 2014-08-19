#ifndef CURVES_TYPED_CURVE_HPP
#define CURVES_TYPED_CURVE_HPP

#include "Curve.hpp"
#include "TypedEvaluator.hpp"

namespace curves {

template<typename CurveConfig>
class TypedCurve : public Curve
{
 public:
  
  /// The value type of the curve.
  typedef typename CurveConfig::ValueType ValueType;

  /// The curve's derivative type.
  typedef typename CurveConfig::DerivativeType DerivativeType;

  typedef TypedEvaluator<CurveConfig> EvaluatorType;
  
  typedef typename EvaluatorType::Ptr EvaluatorTypePtr;

  TypedCurve() { }
  virtual ~TypedCurve() { }

  /// \brief The dimension of the underlying manifold
  virtual size_t dim() const = 0;

  /// \name Methods to evaluate the curve.
  ///@{

  /// Evaluate the ambient space of the curve.
  virtual ValueType evaluateVector(Time time) = 0;
  
  /// Evaluate the curve derivatives.
  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) = 0;

  /// \brief Get an evaluator at this time.
  virtual EvaluatorTypePtr getTypedEvaluator(Time time) = 0;

  ///@}


  /// \name Methods to grow or shrink the curve.
  ///@{

  /// Extend the curve into the past so that it can be evaluated
  /// at these times. The times should be less than 
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values) = 0;

  /// \brief Fit a new curve to these data points
  ///
  /// Underneath the curve should have some default policy for fitting
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values) = 0;

  ///@}
};

} // namespace curves


#endif /* CURVES_TYPED_CURVE_HPP */
