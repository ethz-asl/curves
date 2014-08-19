#ifndef CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP

namespace curves {

class LinearInterpolationVectorSpaceEvaluator {
 public:
  LinearInterpolationVectorSpaceEvaluator();
  virtual ~LinearInterpolationVectorSpaceEvaluator();

  virtual void getKeys(std::vector<Key>& outKeys);
  
  virtual void getCoefficients(std::vector<Key>& outCoefficients);

  /// Evaluate the ambient space of the curve.
  virtual Eigen::VectorXd evaluateVector();
  
  /// Evaluate the curve derivatives.
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder);

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateVector(const std::vector<Coefficient>& coefficients);
  
  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder, const std::vector<Coefficient>& coefficients);

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateVectorAndJacobian(const std::vector<Coefficient>& coefficients, std::vector<Eigen::MatrixXd>& outJacobian);
  
  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivativeAndJacobian(unsigned derivativeOrder, const std::vector<Coefficient>& coefficients, std::vector<Eigen::MatrixXd>& outJacobian);

};

} // namespace curves


#endif /* CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_EVALUATOR_HPP */
