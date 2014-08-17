#ifndef CURVES_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_VECTOR_SPACE_EVALUATOR_HPP

#include "Evaluator.hpp"

namespace curves {

class VectorSpaceEvaluator : public Evaluator
{
 public:
  typedef std::shared_ptr<VectorSpaceEvaluator> Ptr;
  
  VectorSpaceEvaluator();
  virtual ~VectorSpaceEvaluator();

  /// Evaluate the ambient space of the curve.
  virtual Eigen::VectorXd evaluateVector() = 0;
  
  /// Evaluate the curve derivatives.
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder) = 0;

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateVector(const std::vector<Coefficient>& coefficients) = 0;
  
  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder, const std::vector<Coefficient>& coefficients) = 0;

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateVectorAndJacobian(const std::vector<Coefficient>& coefficients, std::vector<Eigen::MatrixXd>& outJacobian) = 0;
  
  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivativeAndJacobian(unsigned derivativeOrder, const std::vector<Coefficient>& coefficients, std::vector<Eigen::MatrixXd>& outJacobian) = 0;

};

} // namespace curves


#endif /* CURVES_VECTOR_SPACE_EVALUATOR_HPP */
