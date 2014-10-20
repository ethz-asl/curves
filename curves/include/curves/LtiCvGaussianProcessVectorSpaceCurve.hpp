#ifndef CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP
#define CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP

#include <curves/GaussianProcessVectorSpaceCurve.hpp>
#include <curves/LtiCvGaussianProcessVectorSpacePrior.hpp>

class LtiCvGaussianProcessVectorSpaceEvaluator; // Forward declaration

namespace curves {

/// \class LtiCvGaussianProcessVectorSpaceCurve
/// \brief A special curve interface for GP curves based on a LTI CV prior.
///
/// A special curve interface for Gaussian-process, vector-space trajectories
/// based on a linear, time-invariant, constant-velocity prior. In order to be
/// interchangable with other typical vector-space curve representations while
/// maintaining sparisty, this interface hides the true underlying ambient curve
/// dimension (position + velocity), in favour of the more typical position-only space.
class LtiCvGaussianProcessVectorSpaceCurve : public GaussianProcessVectorSpaceCurve {
 public:
  /// \brief Parent class
  typedef GaussianProcessVectorSpaceCurve Parent;

  /// \brief The value type of the curve.
  typedef Parent::ValueType ValueType;

  /// \brief The derivative type of the curve.
  typedef Parent::DerivativeType DerivativeType;

  /// \brief The evaluator type of the curve.
  typedef Parent::EvaluatorType EvaluatorType;

  /// \brief The evaluator type pointer.
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;

  /// Constructor to make a Gaussian process curve based on a prior with a
  /// linear, time-invariant, constant-velocity model by specifiying the
  /// stationary power spectral density matrix.
  LtiCvGaussianProcessVectorSpaceCurve(Eigen::MatrixXd stationaryPowerSpectralDensity);
  virtual ~LtiCvGaussianProcessVectorSpaceCurve();

  /// Initialize the linear time invariant, constant velocity prior.
  void initialize(Time initialTime, Eigen::VectorXd initialMean, Eigen::MatrixXd initialInverseCovariance);

  /// Check if the prior has been properly initialized.
  bool isInitialized() const;

  /// Append an exogenous input.
  /// Exogenous inputs create a step function that is assumed to hold constant from last known value.
  void addExogenousInput(Time time, const ValueType& value);

  /// Set the discrete time exogenous inputs.
  /// Exogenous inputs create a step function that is assumed to hold constant from last known value.
  void setExogenousInputs(const std::vector<Time>& times, const std::vector<ValueType>& values);

  /// Evaluate the ambient space of the curve.
  virtual Eigen::VectorXd evaluate(Time time) const;

  /// Evaluate the curve derivatives.
  virtual Eigen::VectorXd evaluateDerivative(Time time, unsigned derivativeOrder) const;

  /// \brief Get an evaluator at this time
  EvaluatorTypePtr getEvaluator(const Time& time) const;

  /// \brief Get the curve dimension.
  virtual size_t dim() const;
};

} // namespace curves


#endif /* CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP */
