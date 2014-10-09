#ifndef CURVES_TYPED_EVALUATOR_HPP
#define CURVES_TYPED_EVALUATOR_HPP

#include "EvaluatorBase.hpp"

namespace curves {

class Coefficients;

template<typename CurveConfig>
class Evaluator : public EvaluatorBase
{
 public:

  /// The value type of the curve.
  typedef typename CurveConfig::ValueType ValueType;

  /// The curve's derivative type.
  typedef typename CurveConfig::DerivativeType DerivativeType;

  /// A pointer type for this evaluator
  typedef boost::shared_ptr< Evaluator<CurveConfig> > Ptr;

  /// A pointer type for this evaluator
  typedef boost::shared_ptr< const Evaluator<CurveConfig> > ConstPtr;

  Evaluator() { }
  virtual ~Evaluator() { }

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluate(const std::vector<Coefficient>& coefficients) const = 0;


  /// Evaluate the ambient space of the curve (functional form) with original coefficients.
  virtual ValueType evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                        std::vector<Eigen::MatrixXd>* jacobians) const = 0;

  /// Evaluate the curve derivatives (functional form).
  virtual ValueType evaluateDerivative(unsigned derivativeOrder, const std::vector<Coefficient>& coefficients) const = 0;
  
  /// Evaluate the curve derivatives (functional form).
  virtual ValueType evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                        const std::vector<Coefficient>& coefficients,
                                                        std::vector<Eigen::MatrixXd>* outJacobian) const = 0;
  /// Evaluate the ambient space of the curve
  virtual ValueType evaluate(const Coefficients& coefficients) const = 0;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluateDerivative(unsigned derivativeOrder,
                                       const Coefficients& coefficients) const = 0;

  /// Get the curve Jacobians.
  /// This is the main interface for GTSAM
  virtual void getJacobians(unsigned derivativeOrder,
                            const Coefficients& coefficients,
                            const Eigen::MatrixXd& chainRule,
                            const std::vector<Eigen::MatrixXd*>& jacobians) const = 0;

  /// Get the maximum derivative order supported by this evaluator.
  size_t getMaximumDerivativeOrder() const;

};

} // namespace 

#endif /* CURVES_TYPED_EVALUATOR_HPP */
