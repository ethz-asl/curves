#ifndef CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP

#include "curves/LtiCvGaussianProcessVectorSpaceCurve.hpp"
#include "curves/GaussianProcessVectorSpaceEvaluator.hpp"

namespace curves {

class Coefficients;

/// \class LtiCvGaussianProcessVectorSpaceEvaluator
/// \brief A special evaluator interface for GP curves based on a LTI CV prior.
///
/// A special evaluator interface for Gaussian-process, vector-space trajectories
/// based on a linear, time-invariant, constant-velocity prior. In order to be
/// interchangable with other typical vector-space evaluators while maintaining
/// sparisty, this interface hides the true underlying ambient curve dimension
/// (position + velocity), in favour of the more typical position-only space.
class LtiCvGaussianProcessVectorSpaceEvaluator : public GaussianProcessVectorSpaceEvaluator {

 public:
  /// \brief Parent class
  typedef GaussianProcessVectorSpaceEvaluator Parent;

  /// \brief The value type of the curve.
  typedef Parent::ValueType ValueType;

  /// \brief The derivative type of the curve.
  typedef Parent::DerivativeType DerivativeType;

  /// Constructor to make a LTI, constant-velocity Gaussian process curve evaluator.
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
};

} // namespace curves


#endif /* CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP */
