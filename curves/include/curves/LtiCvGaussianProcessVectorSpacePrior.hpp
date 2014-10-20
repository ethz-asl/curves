#ifndef CURVES_LINEAR_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_LINEAR_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "LtiSdeGaussianProcessVectorSpacePrior.hpp"

namespace curves {

/// \class LtiCvGaussianProcessVectorSpacePrior
///
/// Interface for vector-space, Gaussian-process-trajectory priors based on
/// linear, time-invariant, stochastic different equation using a constant
/// velocity model.
class LtiCvGaussianProcessVectorSpacePrior : public LtiSdeGaussianProcessVectorSpacePrior {
 public:
  /// \brief Parent class
  typedef LtiSdeGaussianProcessVectorSpacePrior Parent;

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
  /// time-invariant, constant-velocity model by specifiying the
  /// stationary power spectral density matrix.
  LtiCvGaussianProcessVectorSpacePrior(Eigen::MatrixXd stationaryPowerSpectralDensity);
  virtual ~LtiCvGaussianProcessVectorSpacePrior();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// To be implemented by SDE-form-specific class
  /// \todo needs better name
  /// \todo this level could implement numerical integration
  virtual Eigen::VectorXd calculateLiftedExogenousInput(Time time1, Time time2) const;

  /// To be implemented by SDE-form-specific class
  /// \todo this level could implement numerical integration
  virtual Eigen::MatrixXd calculateStateTransitionMatrix(Time time1, Time time2) const;

  /// To be implemented by SDE-form-specific class
  /// \todo needs better name
  /// \todo this level could implement numerical integration
  virtual Eigen::MatrixXd calculateLiftedCovarianceMatrix(Time time1, Time time2) const;
  virtual Eigen::MatrixXd calculateInverseLiftedCovarianceMatrix(Time time1, Time time2) const;
};

} // namespace curves

#endif /* CURVES_LINEAR_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
