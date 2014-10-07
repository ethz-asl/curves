#ifndef CURVES_LINEAR_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_LINEAR_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "LtiSdeGaussianProcessVectorSpacePrior.hpp"

namespace curves {

class LtiCvGaussianProcessVectorSpacePrior : public LtiSdeGaussianProcessVectorSpacePrior {
 public:
  typedef LtiSdeGaussianProcessVectorSpacePrior Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::EvaluatorType EvaluatorType;
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;
  typedef Parent::CurveCoefficientManagerType CurveCoefficientManagerType;

  /// \brief Initialize with the dimension of the vector space
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
