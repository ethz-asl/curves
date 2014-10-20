#include <curves/GaussianProcessVectorSpaceEvaluator.hpp>
#include <curves/Coefficients.hpp>
#include <iostream>

using namespace std;

namespace curves {

GaussianProcessVectorSpaceEvaluator::GaussianProcessVectorSpaceEvaluator(const GaussianProcessVectorSpaceCurve& curve,
                                                                         const Time& time) : prior_(curve.getPrior()), queryTime_(time) {
  std::vector<KeyCoefficientTime*> coeffs;
  curve.appendCoefficientsAt(queryTime_, &coeffs);

  for (std::vector<KeyCoefficientTime*>::const_iterator it = coeffs.begin(); it != coeffs.end(); ++it) {
    keys_.push_back((*it)->key);
    keyTimes_.push_back((*it)->time);
  }
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

std::vector<Key>::const_iterator GaussianProcessVectorSpaceEvaluator::keyBegin() const {
  return keys_.begin();
}

std::vector<Key>::const_iterator GaussianProcessVectorSpaceEvaluator::keyEnd() const {
  return keys_.end();
}

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluate(
    const std::vector<Coefficient>& coefficients) const {
  CHECK_EQ(coefficients.size(), keys_.size()) << "coefficients and keys do not have matching size.";

  // Get interpolation information from prior
  std::vector<Eigen::VectorXd> outMeanAtKeyTimes;
  std::vector<Eigen::MatrixXd> outKtKinvMats;
  Eigen::VectorXd mean = prior_->evaluateAndInterpMatrices(queryTime_, keyTimes_, &outMeanAtKeyTimes, &outKtKinvMats);

  // Calculate the interpolation
  GaussianProcessVectorSpaceEvaluator::ValueType result = mean;
  for (unsigned int i = 0; i < keys_.size(); i++) {
    result = result + outKtKinvMats.at(i) * (coefficients.at(i).getValue() - outMeanAtKeyTimes.at(i));
  }
  return result;
}

GaussianProcessVectorSpaceEvaluator::ValueType GaussianProcessVectorSpaceEvaluator::evaluate(
    const Coefficients& coefficients) const {

  // Get interpolation information from prior
  std::vector<Eigen::VectorXd> outMeanAtKeyTimes;
  std::vector<Eigen::MatrixXd> outKtKinvMats;
  Eigen::VectorXd mean = prior_->evaluateAndInterpMatrices(queryTime_, keyTimes_, &outMeanAtKeyTimes, &outKtKinvMats);

  // Calculate the interpolation
  GaussianProcessVectorSpaceEvaluator::ValueType result = mean;
  for (unsigned int i = 0; i < keys_.size(); i++) {
    result = result + outKtKinvMats.at(i) * (coefficients.get(keys_[i]).getValue() - outMeanAtKeyTimes.at(i));
  }
  return result;
}

// TODO - this function is pretty inefficient! would be better to just use evalAndJacobian for GPs, unless we save Jacs.
void GaussianProcessVectorSpaceEvaluator::getJacobians(unsigned derivativeOrder,
                                                       const Coefficients& /* coefficients */,
                                                       const Eigen::MatrixXd& chainRule,
                                                       const std::vector<Eigen::MatrixXd*>& jacobians) const {
  // TODO(Sean) implement velocity
  CHECK_EQ(derivativeOrder, 0);
  CHECK_EQ(jacobians.size(), keys_.size()) << "number of jacobians and keys do not match.";
  for (unsigned i = 0; i < jacobians.size(); i++) {
    CHECK_NOTNULL(jacobians[i]);
  }

  // Get interpolation information from prior
  std::vector<Eigen::VectorXd> outMeanAtKeyTimes;
  std::vector<Eigen::MatrixXd> outKtKinvMats;
  prior_->evaluateAndInterpMatrices(queryTime_, keyTimes_, &outMeanAtKeyTimes, &outKtKinvMats);

  //TODO check matrix sizes should be chainRule.rows() x coefficient.ndim()
  for (unsigned int i = 0; i < keys_.size(); i++) {
    *(jacobians[i]) += chainRule * outKtKinvMats[i];
  }
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
