#ifndef CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP

#include "Evaluator.hpp"
#include "VectorSpaceConfig.hpp"
#include "GaussianProcessVectorSpaceCurve.hpp"

namespace curves {

class Coefficients;

class GaussianProcessVectorSpaceEvaluator : public Evaluator<VectorSpaceConfig> {

 public:

  typedef Evaluator<VectorSpaceConfig> Parent;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::ValueType ValueType;

  GaussianProcessVectorSpaceEvaluator(const GaussianProcessVectorSpaceCurve& curve, const Time& time);
  virtual ~GaussianProcessVectorSpaceEvaluator();

  virtual void getKeys(std::vector<Key> *outKeys) const;

  void appendKeys(std::vector<Key> *outKeys) const;

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
                            const std::vector<Eigen::MatrixXd*>* jacobians) const;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluateDerivative(unsigned derivativeOrder,
                                       const Coefficients& coefficients) const;

 private:

  //const GaussianProcessVectorSpaceCurve& curve_; // todo: remove?
  std::vector<Key> keys_;
  std::vector<Coefficient> coefficients_;
  std::vector<Coefficient> priorMeanCoefficients_;
  ValueType priorMeanEval_;
  std::vector<Eigen::MatrixXd> interpMatEvals_;
  size_t dimension_;

};

} // namespace curves


#endif /* CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP */
