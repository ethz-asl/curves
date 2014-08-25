#ifndef CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP

#include "Evaluator.hpp"
#include "VectorSpaceConfig.hpp"
#include "LinearInterpolationVectorSpaceCurve.hpp"

namespace curves {

class LinearInterpolationVectorSpaceEvaluator : public Evaluator<VectorSpaceConfig> {
 private:

  const LinearInterpolationVectorSpaceCurve& curve_;
  std::vector<Key> keys_;
  std::vector<Coefficient> coefficients_;
  double alpha_;
  double oneMinusAlpha_;
  size_t dimension_;
  std::vector<Eigen::MatrixXd> jacobians_;

 public:
  typedef Evaluator<VectorSpaceConfig> Parent;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::ValueType ValueType;

  LinearInterpolationVectorSpaceEvaluator(const LinearInterpolationVectorSpaceCurve& curve, const Time& time);
  virtual ~LinearInterpolationVectorSpaceEvaluator();

  virtual void getKeys(std::vector<Key> *outKeys) const;

  virtual void getCoefficients(std::vector<Coefficient>* outCoefficients) const;

  /// Evaluate the ambient space of the curve (functional form) with original coefficients.
  virtual ValueType evaluate() const;

  /// Evaluate the ambient space of the curve (functional form) by specifying new coefficients.
  virtual ValueType evaluate(const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the ambient space of the curve (functional form) with original coefficients.
  virtual ValueType evaluateAndJacobians(std::vector<Eigen::MatrixXd>* outJacobians) const;

  /// Evaluate the ambient space of the curve (functional form) by specifying new coefficients.
  virtual ValueType evaluateAndJacobians(const std::vector<Coefficient>& coefficients,
                                        std::vector<Eigen::MatrixXd>* outJacobians) const;

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
