#ifndef CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP

#include "curves/LtiCvGaussianProcessVectorSpaceCurve.hpp"
#include "curves/GaussianProcessVectorSpaceEvaluator.hpp"

namespace curves {

class Coefficients;

class LtiCvGaussianProcessVectorSpaceEvaluator : public GaussianProcessVectorSpaceEvaluator {

 public:

  typedef GaussianProcessVectorSpaceEvaluator Parent;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::ValueType ValueType;

  LtiCvGaussianProcessVectorSpaceEvaluator(const LtiCvGaussianProcessVectorSpaceCurve& curve, const Time& time);
  virtual ~LtiCvGaussianProcessVectorSpaceEvaluator();

  /// Evaluate the ambient space of the curve (functional form) by specifying new coefficients.
  virtual ValueType evaluate(const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder,
                                             const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                                    std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                        const std::vector<Coefficient>& coefficients,
                                                        std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluate(const Coefficients& coefficients) const;

  /// Get the curve Jacobians.
  /// This is the main interface for GTSAM
  virtual void getJacobians(unsigned derivativeOrder,
                            const Coefficients& coefficients,
                            const Eigen::MatrixXd& chainRule,
                            const std::vector<Eigen::MatrixXd*>& jacobians) const;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluateDerivative(unsigned derivativeOrder,
                                       const Coefficients& coefficients) const;

 private:


};

} // namespace curves


#endif /* CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP */
