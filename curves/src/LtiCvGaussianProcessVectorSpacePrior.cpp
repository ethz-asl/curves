#include <curves/LtiCvGaussianProcessVectorSpacePrior.hpp>
#include <iostream>

namespace curves {

LtiCvGaussianProcessVectorSpacePrior::LtiCvGaussianProcessVectorSpacePrior(Eigen::MatrixXd stationaryPowerSpectralDensity) :
    LtiSdeGaussianProcessVectorSpacePrior((Eigen::MatrixXd(2*stationaryPowerSpectralDensity.rows(),2*stationaryPowerSpectralDensity.rows()) <<
                                            Eigen::MatrixXd::Zero(stationaryPowerSpectralDensity.rows(),stationaryPowerSpectralDensity.rows()),
                                            Eigen::MatrixXd::Identity(stationaryPowerSpectralDensity.rows(),stationaryPowerSpectralDensity.rows()),
                                            Eigen::MatrixXd::Zero(stationaryPowerSpectralDensity.rows(),stationaryPowerSpectralDensity.rows()),
                                            Eigen::MatrixXd::Zero(stationaryPowerSpectralDensity.rows(),stationaryPowerSpectralDensity.rows())).finished(),
                                          (Eigen::MatrixXd(2*stationaryPowerSpectralDensity.rows(),stationaryPowerSpectralDensity.rows()) <<
                                            Eigen::MatrixXd::Zero(stationaryPowerSpectralDensity.rows(),stationaryPowerSpectralDensity.rows()),
                                            Eigen::MatrixXd::Identity(stationaryPowerSpectralDensity.rows(),stationaryPowerSpectralDensity.rows())).finished(),
                                          stationaryPowerSpectralDensity) {}

LtiCvGaussianProcessVectorSpacePrior::~LtiCvGaussianProcessVectorSpacePrior() {}

void LtiCvGaussianProcessVectorSpacePrior::print(const std::string& str) const {
  std::cout << "====================================" << std::endl;
  std::cout << "======= LIN CV VEL GP PRIOR ========" << std::endl;
  std::cout << str << std::endl;
  LtiSdeGaussianProcessVectorSpacePrior::print("=> EXTENDED BY LIN CV VEL GP PRIOR");
}

/// To be implemented by SDE-form-specific class
/// \todo needs better name
Eigen::VectorXd LtiCvGaussianProcessVectorSpacePrior::calculateLiftedExogenousInput(Time time1, Time time2) const {
  CHECK_GT(time2, time1) << "This is probably a mistake, integrating the prior backwards in time is unusual.";
  return Eigen::VectorXd::Zero(this->dim());
}

/// To be implemented by SDE-form-specific class
Eigen::MatrixXd LtiCvGaussianProcessVectorSpacePrior::calculateStateTransitionMatrix(Time time1, Time time2) const {
  unsigned halfDim = this->dim()/2;
  return (Eigen::MatrixXd(this->dim(),this->dim()) << Eigen::MatrixXd::Identity(halfDim,halfDim), Eigen::MatrixXd::Identity(halfDim,halfDim)*(time1-time2),
                                                      Eigen::MatrixXd::Zero(halfDim,halfDim),     Eigen::MatrixXd::Identity(halfDim,halfDim)).finished();
}

/// To be implemented by SDE-form-specific class
/// \todo needs better name
Eigen::MatrixXd LtiCvGaussianProcessVectorSpacePrior::calculateLiftedCovarianceMatrix(Time time1, Time time2) const {
  CHECK_GT(time2, time1) << "This is probably a mistake, integrating the prior backwards in time is unusual.";
  double delt = time2 - time1;
  double delt2 = delt*delt;
  const Eigen::MatrixXd& Qc = getPowerSpectralDensityMatrix();
  return (Eigen::MatrixXd(this->dim(),this->dim()) << Qc*(1.0/3.0)*delt*delt2, Qc*(1.0/2.0)*delt2,
                                                      Qc*(1.0/2.0)*delt2,      Qc*delt).finished();
}
Eigen::MatrixXd LtiCvGaussianProcessVectorSpacePrior::calculateInverseLiftedCovarianceMatrix(Time time1, Time time2) const {
  CHECK_GT(time2, time1) << "This is probably a mistake, integrating the prior backwards in time is unusual.";
  double invdelt = 1.0/(time2 - time1);
  double invdelt2 = invdelt*invdelt;
  const Eigen::MatrixXd& invQc = getInversePowerSpectralDensityMatrix();
  return (Eigen::MatrixXd(this->dim(),this->dim()) << invQc*12*invdelt*invdelt2, invQc*-6*invdelt2,
                                                      invQc*-6*invdelt2,         invQc*4*invdelt).finished();
}

} // namespace curves
