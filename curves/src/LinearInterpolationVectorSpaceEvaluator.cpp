#include <curves/LinearInterpolationVectorSpaceEvaluator.hpp>
#include <iostream>
using namespace std;

namespace curves {

LinearInterpolationVectorSpaceEvaluator::LinearInterpolationVectorSpaceEvaluator(const LinearInterpolationVectorSpaceCurve& curve,
                                                                                 const Time& time) : curve_(curve) {
  KeyCoefficientTime *coeff0, *coeff1;
  curve.getCoefficientsAt(time, &coeff0, &coeff1);
  keys_.push_back(coeff0->key);
  keys_.push_back(coeff1->key);
  coefficients_.push_back(coeff0->coefficient);
  coefficients_.push_back(coeff1->coefficient);

  // Compute alpha_ and oneMinusAlpha_
  Time dt = coeff1->time - coeff0->time;
  CHECK_NE(dt,0) << "requested division by 0";
  Time t = time - coeff0->time;
  alpha_ = double(t)/double(dt);
  oneMinusAlpha_ = 1 - alpha_;

  dimension_ = curve.dim();
}

LinearInterpolationVectorSpaceEvaluator::~LinearInterpolationVectorSpaceEvaluator() {}

void LinearInterpolationVectorSpaceEvaluator::getKeys(std::vector<Key> *outKeys) const {
  CHECK_NOTNULL(outKeys);
  *outKeys = keys_;
}

void LinearInterpolationVectorSpaceEvaluator::appendKeys(std::vector<Key> *outKeys) const {
  CHECK_NOTNULL(outKeys);
  outKeys->insert(outKeys->end(), keys_.begin(), keys_.end());
}

void LinearInterpolationVectorSpaceEvaluator::getCoefficients(std::vector<Coefficient> *outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  *outCoefficients = coefficients_;
}

void LinearInterpolationVectorSpaceEvaluator::appendCoefficients(std::vector<Coefficient> *outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  outCoefficients->insert(outCoefficients->end(), coefficients_.begin(), coefficients_.end());
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluate() const {
  return coefficients_[0].getValue() * oneMinusAlpha_ + coefficients_[1].getValue() * alpha_;
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluate(
    const std::vector<Coefficient>& coefficients) const {

  return coefficients[0].getValue() * oneMinusAlpha_ + coefficients[1].getValue() * alpha_;
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluateAndJacobians(
    const Eigen::MatrixXd& chainRule,
    const std::vector<Eigen::MatrixXd*>& jacobians) const {

  CHECK_NOTNULL(jacobians[0]);
  CHECK_NOTNULL(jacobians[1]);
  *(jacobians[0]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * oneMinusAlpha_;
  *(jacobians[1]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * alpha_;
  return evaluate();
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluateAndJacobians(
    const std::vector<Coefficient>& coefficients,
    const Eigen::MatrixXd& chainRule,
    const std::vector<Eigen::MatrixXd*>& jacobians) const {

  CHECK_NOTNULL(jacobians[0]);
  CHECK_NOTNULL(jacobians[1]);
  *(jacobians[0]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * oneMinusAlpha_;
  *(jacobians[1]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * alpha_;
  return evaluate(coefficients);
}

/// Evaluate the curve derivatives.
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder) const {
  // \todo Abel and Renaud
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
