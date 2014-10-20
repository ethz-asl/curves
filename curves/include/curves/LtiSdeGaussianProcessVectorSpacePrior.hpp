#ifndef CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "LinearSdeGaussianProcessVectorSpacePrior.hpp"

namespace curves {

/// \class LtiSdeGaussianProcessVectorSpacePrior
///
/// Interface for vector-space, Gaussian-process-trajectory priors based on
/// linear, time-invariant, stochastic different equation using a constant
/// velocity model.
class LtiSdeGaussianProcessVectorSpacePrior : public LinearSdeGaussianProcessVectorSpacePrior {
 public:
  /// \brief Parent class
  typedef LinearSdeGaussianProcessVectorSpacePrior Parent;

  /// \brief The value type of the curve.
  typedef Parent::ValueType ValueType;

  /// \brief The derivative type of the curve.
  typedef Parent::DerivativeType DerivativeType;

  /// \brief The evaluator type of the curve.
  typedef Parent::EvaluatorType EvaluatorType;

  /// \brief The evaluator type pointer.
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;

  /// \brief The coefficient manager type the GP curve should use.
  typedef Parent::CurveCoefficientManagerType CurveCoefficientManagerType;

  /// Constructor to make a Gaussian process prior based on a linear,
  /// time-invariant, stochastic differential equation by specifiying the
  /// drift and diffusion coefficients, and the power spectral density.
  LtiSdeGaussianProcessVectorSpacePrior(Eigen::MatrixXd ltiDriftMatrix, Eigen::MatrixXd ltiDiffusionMatrix, Eigen::MatrixXd stationaryPowerSpectralDensity);
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
  virtual Eigen::MatrixXd calculateInverseLiftedCovarianceMatrix(Time time1, Time time2) const = 0;

  /// Get the drift coefficient matrix
  const Eigen::MatrixXd& getDriftMatrix() const {return ltiDriftMatrix_;}

  /// Get the diffusion coefficient matrix
  const Eigen::MatrixXd& getDiffusionMatrix() const {return ltiDiffusionMatrix_;}

 private:

  /// Linear time-invariant drift matrix
  ///   --- see F in equation X, Anderson et al. (TBD)
  const Eigen::MatrixXd ltiDriftMatrix_;

  /// Linear time-invariant diffusion matrix
  ///   --- see L in equation X, Anderson et al. (TBD)
  const Eigen::MatrixXd ltiDiffusionMatrix_;
};

} // namespace curves

#endif /* CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
