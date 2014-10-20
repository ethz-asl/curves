#include <curves/LtiSdeGaussianProcessVectorSpacePrior.hpp>
#include <iostream>

namespace curves {

LtiSdeGaussianProcessVectorSpacePrior::LtiSdeGaussianProcessVectorSpacePrior(Eigen::MatrixXd ltiDriftMatrix, Eigen::MatrixXd ltiDiffusionMatrix, Eigen::MatrixXd stationaryPowerSpectralDensity) :
    LinearSdeGaussianProcessVectorSpacePrior(ltiDriftMatrix.rows(), stationaryPowerSpectralDensity), ltiDriftMatrix_(ltiDriftMatrix), ltiDiffusionMatrix_(ltiDiffusionMatrix) {
  CHECK_EQ(ltiDriftMatrix_.rows(), ltiDriftMatrix_.cols()) << "Linear time-invariant drift matrix not square.";
  CHECK_EQ(ltiDriftMatrix_.rows(), ltiDiffusionMatrix_.rows()) << "Drift matrix and diffusion matrix do not have the same `output' dimension.";
  CHECK_EQ(ltiDiffusionMatrix_.cols(), stationaryPowerSpectralDensity.rows()) << "Diffusion matrix columns and stationary power-spectral density matrix rows are not equal.";
}

LtiSdeGaussianProcessVectorSpacePrior::~LtiSdeGaussianProcessVectorSpacePrior() {}

void LtiSdeGaussianProcessVectorSpacePrior::print(const std::string& str) const {
  std::cout << "=================================" << std::endl;
  std::cout << "======= LTI SDE GP PRIOR ========" << std::endl;
  std::cout << str << std::endl;
  std::cout << "static drift matrix: " << std::endl << ltiDriftMatrix_ << std::endl;
  std::cout << "static diffusion matrix: " << std::endl << ltiDiffusionMatrix_ << std::endl;
  LinearSdeGaussianProcessVectorSpacePrior::print("=> EXTENDED BY LTI SDE GP PRIOR");
}

} // namespace curves
