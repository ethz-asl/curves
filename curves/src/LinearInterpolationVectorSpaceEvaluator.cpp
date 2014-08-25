#include <curves/LinearInterpolationVectorSpaceEvaluator.hpp>
#include <iostream>
using namespace std;

namespace curves {

LinearInterpolationVectorSpaceEvaluator::LinearInterpolationVectorSpaceEvaluator(const LinearInterpolationVectorSpaceCurve& curve,
                                                                                 const Time& time) : curve_(curve) {
  vector<KeyCoefficientTime*> coefficients;
  KeyCoefficientTime coeff0, coeff1;
  coefficients.push_back(&coeff0);
  coefficients.push_back(&coeff1);
  curve.getCoefficientsAt(time, coefficients);
  keys_.push_back(coefficients[0]->key);
  keys_.push_back(coefficients[1]->key);
  coefficients_.push_back(coefficients[0]->coefficient);
  coefficients_.push_back(coefficients[1]->coefficient);

  // Compute alpha_ and oneMinusAlpha_
  Time dt = coefficients[1]->time - coefficients[0]->time;
  Time t = time - coefficients[0]->time;
  alpha_ = double(t)/double(dt);
  oneMinusAlpha_ = 1 - alpha_;

  dimension_ = curve.dim();

  // Compute the Jacobians
  // In the linear case, the jacobians are independent of the coefficients value.
  // J0 and J1 are diagonal matrices with elements equal to oneMinusAlpha_ and alpha_ respectively.
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dimension_,dimension_);
  jacobians_.push_back(I*oneMinusAlpha_);
  jacobians_.push_back(I*alpha_);
}

LinearInterpolationVectorSpaceEvaluator::~LinearInterpolationVectorSpaceEvaluator() {}

void LinearInterpolationVectorSpaceEvaluator::getKeys(std::vector<Key> *outKeys) const {
  *outKeys = keys_;
}

void LinearInterpolationVectorSpaceEvaluator::getCoefficients(std::vector<Coefficient> *outCoefficients) const {
  *outCoefficients = coefficients_;
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluate() const {
  return coefficients_[0].getValue() * oneMinusAlpha_ + coefficients_[1].getValue() * alpha_;
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluate(const std::vector<Coefficient>& coefficients) const {
  return coefficients[0].getValue() * oneMinusAlpha_ + coefficients[1].getValue() * alpha_;
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluateAndJacobians(std::vector<Eigen::MatrixXd>* outJacobians) const {
  CHECK_NOTNULL(outJacobians);
  *outJacobians = jacobians_;
  return evaluate();
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluateAndJacobians(const std::vector<Coefficient>& coefficients,
                                                                                                                 std::vector<Eigen::MatrixXd>* outJacobians) const {
  CHECK_NOTNULL(outJacobians);
  *outJacobians = jacobians_;
  return evaluate(coefficients);
}

/// Evaluate the ambient space of the curve.
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateVector() const {
  // \todo Abel and Renaud -- Needed?
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives.
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the ambient space of the curve (functional form).
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateVector(const std::vector<Coefficient>& coefficients) const {
  // \todo Abel and Renaud -- Needed?
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives (functional form).
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                            const std::vector<Coefficient>& coefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the ambient space of the curve (functional form).
LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluateVectorAndJacobian(const std::vector<Coefficient>& coefficients,
                                                                                                                      std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives (functional form).
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                                                       const std::vector<Coefficient>& coefficients,
                                                                                       std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

} // namespace curves
