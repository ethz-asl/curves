#ifndef CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "LinearSdeGaussianProcessVectorSpacePrior.hpp"
#include "HermiteCoefficientManager.hpp"

namespace curves {

class LtiSdeGaussianProcessVectorSpacePrior : public LinearSdeGaussianProcessVectorSpacePrior {
 public:
  typedef LinearSdeGaussianProcessVectorSpacePrior Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::EvaluatorType EvaluatorType;
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;
  typedef Parent::CurveCoefficientManagerType CurveCoefficientManagerType;

  /// \brief Initialize with the dimension of the vector space
  LtiSdeGaussianProcessVectorSpacePrior(Eigen::MatrixXd ltiDriftMatrix_, Eigen::MatrixXd ltiDiffusionMatrix_, Eigen::MatrixXd stationaryPowerSpectralDensity);
  virtual ~LtiSdeGaussianProcessVectorSpacePrior();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// To be implemented by SDE-form-specific class
  /// \todo needs better name
  /// \todo this level could implement numerical integration
  virtual Eigen::VectorXd calculateLiftedExogenousInput(Time time1, Time time2) const = 0;

  /// To be implemented by SDE-form-specific class
  /// \todo this level could implement numerical integration
  virtual Eigen::MatrixXd calculateStateTransitionMatrix(Time time1, Time time2) const = 0;

  /// To be implemented by SDE-form-specific class
  /// \todo needs better name
  /// \todo this level could implement numerical integration
  virtual Eigen::MatrixXd calculateLiftedCovarianceMatrix(Time time1, Time time2) const = 0;

 private:

  /// Linear time-invariant drift matrix
  ///   --- see F in equation X, Anderson et al. (TBD)
  const Eigen::MatrixXd ltiDriftMatrix__;

  /// Linear time-invariant diffusion matrix
  ///   --- see L in equation X, Anderson et al. (TBD)
  const Eigen::MatrixXd ltiDiffusionMatrix__;

};

} // namespace curves

#endif /* CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
