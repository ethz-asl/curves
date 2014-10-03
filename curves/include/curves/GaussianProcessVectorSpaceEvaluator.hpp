#ifndef CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP

#include "Evaluator.hpp"
#include "VectorSpaceConfig.hpp"
#include "GaussianProcessVectorSpaceCurve.hpp"

namespace curves {

class GaussianProcessVectorSpaceEvaluator : public Evaluator<VectorSpaceConfig> {

 public:

  typedef Evaluator<VectorSpaceConfig> Parent;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::ValueType ValueType;

  GaussianProcessVectorSpaceEvaluator(const GaussianProcessVectorSpaceCurve& curve, const Time& time);
  virtual ~GaussianProcessVectorSpaceEvaluator();

  virtual void getKeys(std::vector<Key> *outKeys) const;

  void appendKeys(std::vector<Key> *outKeys) const;

  virtual void getCoefficients(std::vector<Coefficient>* outCoefficients) const;

  void appendCoefficients(std::vector<Coefficient> *outCoefficients) const;

  /// Evaluate the ambient space of the curve (functional form) with original coefficients.
  virtual ValueType evaluate() const;

  /// Evaluate the ambient space of the curve (functional form) by specifying new coefficients.
  virtual ValueType evaluate(const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the ambient space of the curve (functional form) with original coefficients.
  virtual ValueType evaluateAndJacobians(const Eigen::MatrixXd& chainRule,
                                         const std::vector<Eigen::MatrixXd*>& jacobians) const;

  virtual ValueType evaluateAndJacobians(const std::vector<Coefficient>& coefficients,
                                         const Eigen::MatrixXd& chainRule,
                                         const std::vector<Eigen::MatrixXd*>& jacobians) const;

  /// Evaluate the curve derivatives.
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder) const;

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

 private:

  const GaussianProcessVectorSpaceCurve& curve_; // todo: remove?
  std::vector<Key> keys_;
  std::vector<Coefficient> coefficients_;
  std::vector<Coefficient> priorMeanCoefficients_;
  ValueType priorMeanEval_;
  std::vector<Eigen::MatrixXd> interpMatEvals_;
  //double alpha_;
  //double oneMinusAlpha_;
  size_t dimension_;

};

} // namespace curves


#endif /* CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP */
