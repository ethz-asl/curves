#ifndef CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP
#define CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP

#include <curves/GaussianProcessVectorSpaceCurve.hpp>
#include <curves/LtiCvGaussianProcessVectorSpacePrior.hpp>

class LtiCvGaussianProcessVectorSpaceEvaluator; // Forward declaration

namespace curves {

class LtiCvGaussianProcessVectorSpaceCurve : public GaussianProcessVectorSpaceCurve {
 public:
  typedef GaussianProcessVectorSpaceCurve::ValueType ValueType;
  typedef GaussianProcessVectorSpaceCurve::DerivativeType DerivativeType;
  typedef GaussianProcessVectorSpaceCurve::EvaluatorType EvaluatorType;
  typedef GaussianProcessVectorSpaceCurve::EvaluatorTypePtr EvaluatorTypePtr;

  /// \brief Initialize with the dimension of the vector space
  LtiCvGaussianProcessVectorSpaceCurve(Eigen::MatrixXd stationaryPowerSpectralDensity);
  virtual ~LtiCvGaussianProcessVectorSpaceCurve();

  void initialize(Time initialTime, Eigen::VectorXd initialMean, Eigen::MatrixXd initialInverseCovariance);
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

  virtual size_t dim() const;

 private:


};

} // namespace curves


#endif /* CURVES_LTI_CONSTANT_VELOCITY_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP */
