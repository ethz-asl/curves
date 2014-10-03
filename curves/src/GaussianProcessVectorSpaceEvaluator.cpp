#include <curves/GaussianProcessVectorSpaceEvaluator.hpp>
#include <iostream>
using namespace std;

namespace curves {

GaussianProcessVectorSpaceEvaluator::GaussianProcessVectorSpaceEvaluator(const GaussianProcessVectorSpaceCurve& curve,
                                                                         const Time& time) : curve_(curve) {
  KeyCoefficientTime *coeff0, *coeff1;
  //curve.getCoefficientsAt(time, &coeff0, &coeff1);
  keys_.push_back(coeff0->key);
  keys_.push_back(coeff1->key);
  coefficients_.push_back(coeff0->coefficient);
  coefficients_.push_back(coeff1->coefficient);

  // for each key, get the coefficients?
  //std::vector<Coefficient> priorMeanCoefficients_;
  //ValueType priorMeanEval_;
  //std::vector<Eigen::MatrixXd> interpMatEvals_;

  // Compute alpha_ and oneMinusAlpha_
  Time dt = coeff1->time - coeff0->time;
  CHECK_NE(dt,0) << "requested division by 0";
  Time t = time - coeff0->time;
  //alpha_ = double(t)/double(dt);
  //oneMinusAlpha_ = 1 - alpha_;

  dimension_ = curve.dim();
}

GaussianProcessVectorSpaceEvaluator::~GaussianProcessVectorSpaceEvaluator() {}

void GaussianProcessVectorSpaceEvaluator::getKeys(std::vector<Key> *outKeys) const {
  CHECK_NOTNULL(outKeys);
  *outKeys = keys_;
}

void GaussianProcessVectorSpaceEvaluator::appendKeys(std::vector<Key> *outKeys) const {
  CHECK_NOTNULL(outKeys);
  outKeys->insert(outKeys->end(), keys_.begin(), keys_.end());
}

void GaussianProcessVectorSpaceEvaluator::getCoefficients(std::vector<Coefficient> *outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  *outCoefficients = coefficients_;
}

void GaussianProcessVectorSpaceEvaluator::appendCoefficients(std::vector<Coefficient> *outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  outCoefficients->insert(outCoefficients->end(), coefficients_.begin(), coefficients_.end());
}

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluate() const {
  //return coefficients_[0].getValue() * oneMinusAlpha_ + coefficients_[1].getValue() * alpha_;
  return coefficients_[0].getValue() + coefficients_[1].getValue();
}

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluate(
    const std::vector<Coefficient>& coefficients) const {
  //return coefficients[0].getValue() * oneMinusAlpha_ + coefficients[1].getValue() * alpha_;
  return coefficients[0].getValue() + coefficients[1].getValue();
}

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluateAndJacobians(
    const Eigen::MatrixXd& chainRule,
    const std::vector<Eigen::MatrixXd*>& jacobians) const {

  CHECK_NOTNULL(jacobians[0]);
  CHECK_NOTNULL(jacobians[1]);
  *(jacobians[0]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * 0.5;
  *(jacobians[1]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * 0.5;
  return evaluate();
}

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluateAndJacobians(
    const std::vector<Coefficient>& coefficients,
    const Eigen::MatrixXd& chainRule,
    const std::vector<Eigen::MatrixXd*>& jacobians) const {

  CHECK_NOTNULL(jacobians[0]);
  CHECK_NOTNULL(jacobians[1]);
  *(jacobians[0]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * 0.5;
  *(jacobians[1]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * 0.5;
  return evaluate(coefficients);
}

/// Evaluate the curve derivatives.

Eigen::VectorXd GaussianProcessVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

/// Evaluate the curve derivatives (functional form).

Eigen::VectorXd GaussianProcessVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                            const std::vector<Coefficient>& coefficients) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

/// Evaluate the ambient space of the curve (functional form).

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluateVectorAndJacobian(const std::vector<Coefficient>& coefficients,
                                                                                                              std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

/// Evaluate the curve derivatives (functional form).

Eigen::VectorXd GaussianProcessVectorSpaceEvaluator::evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                                                   const std::vector<Coefficient>& coefficients,
                                                                                   std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

} // namespace curves
