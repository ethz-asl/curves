#include <curves/LtiCvGaussianProcessVectorSpaceCurve.hpp>
#include <curves/LtiCvGaussianProcessVectorSpaceEvaluator.hpp>

namespace curves {

LtiCvGaussianProcessVectorSpaceCurve::LtiCvGaussianProcessVectorSpaceCurve(Eigen::MatrixXd stationaryPowerSpectralDensity) :
  GaussianProcessVectorSpaceCurve(boost::shared_ptr<GaussianProcessVectorSpacePrior>( new LtiCvGaussianProcessVectorSpacePrior(stationaryPowerSpectralDensity) )) {}

LtiCvGaussianProcessVectorSpaceCurve::~LtiCvGaussianProcessVectorSpaceCurve() {}

void LtiCvGaussianProcessVectorSpaceCurve::initialize(Time initialTime, Eigen::VectorXd initialMean, Eigen::MatrixXd initialInverseCovariance) {
  boost::static_pointer_cast<LtiCvGaussianProcessVectorSpacePrior>(this->getPrior())->initialize(initialTime, initialMean, initialInverseCovariance);
}

bool LtiCvGaussianProcessVectorSpaceCurve::isInitialized() const {
  return boost::static_pointer_cast<LtiCvGaussianProcessVectorSpacePrior>(this->getPrior())->isInitialized();
}

/// Append an exogenous input.
/// Exogenous inputs create a step function that is assumed to hold constant from last known value.
void LtiCvGaussianProcessVectorSpaceCurve::addExogenousInput(Time time, const ValueType& value) {
  return boost::static_pointer_cast<LtiCvGaussianProcessVectorSpacePrior>(this->getPrior())->addExogenousInput(time, value);
}

/// Set the discrete time exogenous inputs.
/// Exogenous inputs create a step function that is assumed to hold constant from last known value.
void LtiCvGaussianProcessVectorSpaceCurve::setExogenousInputs(const std::vector<Time>& times, const std::vector<ValueType>& values) {
  return boost::static_pointer_cast<LtiCvGaussianProcessVectorSpacePrior>(this->getPrior())->setExogenousInputs(times, values);
}

Eigen::VectorXd LtiCvGaussianProcessVectorSpaceCurve::evaluate(Time time) const {
  return (Eigen::MatrixXd(this->getPrior()->dim()/2,this->getPrior()->dim()) <<
           Eigen::MatrixXd::Identity(this->getPrior()->dim()/2,this->getPrior()->dim()/2),
           Eigen::MatrixXd::Zero(this->getPrior()->dim()/2,this->getPrior()->dim()/2)).finished() *
           this->GaussianProcessVectorSpaceCurve::evaluate(time);
}

Eigen::VectorXd LtiCvGaussianProcessVectorSpaceCurve::evaluateDerivative(Time time, unsigned derivativeOrder) const {
  // \todo Sean
  if (derivativeOrder == 0) {
    return this->evaluate(time);
  } else {
    return (Eigen::MatrixXd(this->getPrior()->dim()/2,this->getPrior()->dim()) <<
            Eigen::MatrixXd::Zero(this->getPrior()->dim()/2,this->getPrior()->dim()/2),
            Eigen::MatrixXd::Identity(this->getPrior()->dim()/2,this->getPrior()->dim()/2)).finished() *
            this->GaussianProcessVectorSpaceCurve::evaluateDerivative(time, derivativeOrder-1);
  }
}

/// \brief Get an evaluator at this time
typename LtiCvGaussianProcessVectorSpaceCurve::EvaluatorTypePtr LtiCvGaussianProcessVectorSpaceCurve::getEvaluator(const Time& time) const {
  boost::shared_ptr< Evaluator<VectorSpaceConfig> > rval( new LtiCvGaussianProcessVectorSpaceEvaluator((*this), time) );
  return rval;
}

size_t LtiCvGaussianProcessVectorSpaceCurve::dim() const {
  return this->GaussianProcessVectorSpaceCurve::dim()/2;
}

} // namespace curves
