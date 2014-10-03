#include <curves/GaussianProcessVectorSpaceEvaluator.hpp>
#include <curves/Coefficients.hpp>
#include <iostream>
using namespace std;

namespace curves {

GaussianProcessVectorSpaceEvaluator::GaussianProcessVectorSpaceEvaluator(const GaussianProcessVectorSpaceCurve& curve,
                                                                         const Time& time) {
  KeyCoefficientTime *coeff0, *coeff1;
  //curve.getCoefficientsAt(time, &coeff0, &coeff1);
  keys_.push_back(coeff0->key);
  keys_.push_back(coeff1->key);

  // for each key, get the coefficients?
  //std::vector<Coefficient> priorMeanCoefficients_;
  //ValueType priorMeanEval_;
  //std::vector<Eigen::MatrixXd> interpMatEvals_;

  // Compute alpha_
  Time dt = coeff1->time - coeff0->time;
  if(dt == 0) {
    //alpha_ = 1.0;
  } else {
    Time t = time - coeff0->time;
    //alpha_ = double(t)/double(dt);
  }

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

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluate(
    const std::vector<Coefficient>& coefficients) const {

  return coefficients[0].getValue() * 0.5 + coefficients[1].getValue() * 0.5;
}

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluate(
    const Coefficients& coefficients) const {

  return coefficients.get(keys_[0]).getValue() * 0.5 + coefficients.get(keys_[1]).getValue() * 0.5;
}

void GaussianProcessVectorSpaceEvaluator::getJacobians(unsigned derivativeOrder,
                                                       const Coefficients& /* coefficients */,
                                                       const Eigen::MatrixXd& chainRule,
                                                       const std::vector<Eigen::MatrixXd*>* jacobians) const {
  CHECK_EQ(derivativeOrder, 0);
  CHECK_NOTNULL(jacobians);
  CHECK_NOTNULL((*jacobians)[0]);
  CHECK_NOTNULL((*jacobians)[1]);
  *((*jacobians)[0]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * 0.5;
  *((*jacobians)[1]) += chainRule * Eigen::MatrixXd::Identity(dimension_,dimension_) * 0.5;
}

///

/// Evaluate the ambient space of the curve
Eigen::VectorXd GaussianProcessVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                        const Coefficients& coefficients) const {
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
Eigen::VectorXd GaussianProcessVectorSpaceEvaluator::evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
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
