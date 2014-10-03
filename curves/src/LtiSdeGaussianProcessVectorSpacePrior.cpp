#include <curves/LtiSdeGaussianProcessVectorSpacePrior.hpp>
#include <iostream>

namespace curves {

LtiSdeGaussianProcessVectorSpacePrior::LtiSdeGaussianProcessVectorSpacePrior(Eigen::MatrixXd ltiDriftMatrix_, Eigen::MatrixXd ltiDiffusionMatrix_, Eigen::MatrixXd stationaryPowerSpectralDensity) :
    LinearSdeGaussianProcessVectorSpacePrior(ltiDriftMatrix_.rows(), stationaryPowerSpectralDensity), ltiDriftMatrix__(ltiDriftMatrix_), ltiDiffusionMatrix__(ltiDiffusionMatrix_) {
  CHECK_EQ(ltiDriftMatrix_.rows(), ltiDriftMatrix_.cols()) << "Linear time-invariant drift matrix not square.";
  CHECK_EQ(ltiDriftMatrix_.rows(), ltiDiffusionMatrix_.rows()) << "Drift matrix and diffusion matrix do not have the same `output' dimension.";
  CHECK_EQ(ltiDiffusionMatrix_.cols(), stationaryPowerSpectralDensity.rows()) << "Diffusion matrix columns and stationary power-spectral density matrix rows are not equal.";
}

LtiSdeGaussianProcessVectorSpacePrior::~LtiSdeGaussianProcessVectorSpacePrior() {}

void LtiSdeGaussianProcessVectorSpacePrior::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=======LTI GAUSSIAN PROCESS PRIOR========" << std::endl;
  std::cout << "=========================================" <<std::endl;
}



} // namespace curves
