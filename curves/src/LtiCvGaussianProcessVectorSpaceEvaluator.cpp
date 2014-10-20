#include <curves/LtiCvGaussianProcessVectorSpaceEvaluator.hpp>
//#include <curves/Coefficients.hpp>
#include <iostream>

using namespace std;

namespace curves {

LtiCvGaussianProcessVectorSpaceEvaluator::LtiCvGaussianProcessVectorSpaceEvaluator(const LtiCvGaussianProcessVectorSpaceCurve& curve,
                                                                                   const Time& time) : GaussianProcessVectorSpaceEvaluator(curve, time) {}

LtiCvGaussianProcessVectorSpaceEvaluator::~LtiCvGaussianProcessVectorSpaceEvaluator() {}

LtiCvGaussianProcessVectorSpaceEvaluator::ValueType LtiCvGaussianProcessVectorSpaceEvaluator::evaluate(
    const std::vector<Coefficient>& coefficients) const {
  size_t priorDim = prior_->dim();
  return (Eigen::MatrixXd(priorDim/2,priorDim) <<
           Eigen::MatrixXd::Identity(priorDim/2,priorDim/2),
           Eigen::MatrixXd::Zero(priorDim/2,priorDim/2)).finished() *
           this->GaussianProcessVectorSpaceEvaluator::evaluate(coefficients);
}

LtiCvGaussianProcessVectorSpaceEvaluator::ValueType LtiCvGaussianProcessVectorSpaceEvaluator::evaluate(
    const Coefficients& coefficients) const {
  size_t priorDim = prior_->dim();
  return (Eigen::MatrixXd(priorDim/2,priorDim) <<
           Eigen::MatrixXd::Identity(priorDim/2,priorDim/2),
           Eigen::MatrixXd::Zero(priorDim/2,priorDim/2)).finished() *
           this->GaussianProcessVectorSpaceEvaluator::evaluate(coefficients);
}

// TODO - this function is pretty inefficient! would be better to just use evalAndJacobian for GPs, unless we save Jacs.
void LtiCvGaussianProcessVectorSpaceEvaluator::getJacobians(unsigned derivativeOrder,
                                                            const Coefficients& /* coefficients */,
                                                            const Eigen::MatrixXd& chainRule,
                                                            const std::vector<Eigen::MatrixXd*>& jacobians) const {
  // TODO(Sean) implement velocity
  CHECK_EQ(derivativeOrder, 0);
  CHECK_EQ(jacobians.size(), keys_.size()) << "number of jacobians and keys do not match.";
  for (unsigned i = 0; i < jacobians.size(); i++) {
    CHECK_NOTNULL(jacobians[i]);
  }
  size_t priorDim = prior_->dim();

  // Get interpolation information from prior
  std::vector<Eigen::VectorXd> outMeanAtKeyTimes;
  std::vector<Eigen::MatrixXd> outKtKinvMats;
  prior_->evaluateAndInterpMatrices(queryTime_, keyTimes_, &outMeanAtKeyTimes, &outKtKinvMats);

  // Create projection matrix
  Eigen::MatrixXd projection(priorDim/2,priorDim);
  projection << Eigen::MatrixXd::Identity(priorDim/2,priorDim/2),
                Eigen::MatrixXd::Zero(priorDim/2,priorDim/2);

  //TODO check matrix sizes should be chainRule.rows() x coefficient.ndim()
  for (unsigned int i = 0; i < keys_.size(); i++) {
    *(jacobians[i]) += chainRule * projection * outKtKinvMats[i];
  }
}

/// Evaluate the ambient space of the curve
Eigen::VectorXd LtiCvGaussianProcessVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                        const Coefficients& coefficients) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

/// Evaluate the curve derivatives (functional form).
Eigen::VectorXd LtiCvGaussianProcessVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                        const std::vector<Coefficient>& coefficients) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

/// Evaluate the ambient space of the curve (functional form).
Eigen::VectorXd LtiCvGaussianProcessVectorSpaceEvaluator::evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                                                          std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

/// Evaluate the curve derivatives (functional form).
Eigen::VectorXd LtiCvGaussianProcessVectorSpaceEvaluator::evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                                                   const std::vector<Coefficient>& coefficients,
                                                                                   std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Sean
  CHECK(false) << "To be implemented by Sean";
}

} // namespace curves
