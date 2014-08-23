#ifndef CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP

#include "Evaluator.hpp"
#include "VectorSpaceConfig.hpp"

namespace curves {

class LinearInterpolationVectorSpaceEvaluator : public Evaluator<VectorSpaceConfig> {

 public:
  typedef Evaluator<VectorSpaceConfig> Parent;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::ValueType ValueType;

  LinearInterpolationVectorSpaceEvaluator() {}
  virtual ~LinearInterpolationVectorSpaceEvaluator() {}

  virtual void getKeys(std::vector<Key> *outKeys) const;

  virtual void getCoefficients(std::vector<Coefficient>* outCoefficients) const;

  /// Evaluate the ambient space of the curve (functional form).
  virtual ValueType evaluate(const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the ambient space of the curve (functional form).
  virtual ValueType evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                        std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the ambient space of the curve.
  virtual Eigen::VectorXd evaluateVector() const;

  /// Evaluate the curve derivatives.
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder) const;

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateVector(const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder,
                                             const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateVectorAndJacobian(const std::vector<Coefficient>& coefficients,
                                                    std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                        const std::vector<Coefficient>& coefficients,
                                                        std::vector<Eigen::MatrixXd>* outJacobian) const;

};

} // namespace curves


#endif /* CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP */
