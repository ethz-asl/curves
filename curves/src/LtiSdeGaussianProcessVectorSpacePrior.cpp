#include <curves/LtiSdeGaussianProcessVectorSpacePrior.hpp>
#include <iostream>

namespace curves {

//LtiSdeGaussianProcessVectorSpacePrior::LtiSdeGaussianProcessVectorSpacePrior(size_t dimension) :
//                                            GaussianProcessVectorSpacePrior(dimension) {}

LtiSdeGaussianProcessVectorSpacePrior::LtiSdeGaussianProcessVectorSpacePrior(Eigen::MatrixXd ltiSdeMatrixF, Eigen::MatrixXd ltiSdeMatrixL, Eigen::MatrixXd ltiSdeMatrixQc) :
    GaussianProcessVectorSpacePrior(ltiSdeMatrixF.rows()), ltiSdeMatrixF_(ltiSdeMatrixF), ltiSdeMatrixL_(ltiSdeMatrixL), ltiSdeMatrixQc_(ltiSdeMatrixQc) {
  CHECK_EQ(ltiSdeMatrixF.rows(), ltiSdeMatrixF.cols()) << "F matrix not square.";
  CHECK_EQ(ltiSdeMatrixF.rows(), ltiSdeMatrixL.rows()) << "F matrix and L matrix not of the same output dimension.";
  CHECK_EQ(ltiSdeMatrixL.cols(), ltiSdeMatrixQc.rows()) << "L matrix columns and Qc matrix rows are not equal.";
  CHECK_EQ(ltiSdeMatrixQc.rows(), ltiSdeMatrixQc.cols()) << "Qc matrix is not square.";
}

LtiSdeGaussianProcessVectorSpacePrior::~LtiSdeGaussianProcessVectorSpacePrior() {}

void LtiSdeGaussianProcessVectorSpacePrior::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=======GAUSSIAN PROCESS PRIOR========" << std::endl;
  std::cout << "=========================================" <<std::endl;
}


void LtiSdeGaussianProcessVectorSpacePrior::getCoefficientsAt(const Time& time,
                                                            Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  (*outCoefficients)[rval0->key] = rval0->coefficient;
  (*outCoefficients)[rval1->key] = rval1->coefficient;
}

void LtiSdeGaussianProcessVectorSpacePrior::getCoefficientsInRange(Time startTime,
                                                                 Time endTime, 
                                                                 Coefficient::Map* outCoefficients) const {
  manager_.getCoefficientsInRange(startTime, endTime, outCoefficients);
}

void LtiSdeGaussianProcessVectorSpacePrior::getCoefficients(Coefficient::Map* outCoefficients) const {
  manager_.getCoefficients(outCoefficients);
}

Time LtiSdeGaussianProcessVectorSpacePrior::getMaxTime() const {
  return manager_.getMaxTime();
}

Time LtiSdeGaussianProcessVectorSpacePrior::getMinTime() const {
  return manager_.getMinTime();
}

/// Append an exogenous input.
/// Exogenous inputs create a step function that is assumed to hold constant from last known value.
void LtiSdeGaussianProcessVectorSpacePrior::addExogenousInput(Time time, const ValueType& value) {
  CHECK(false) << "Not implemented, todo Sean";
}

/// Set the discrete time exogenous inputs.
/// Exogenous inputs create a step function that is assumed to hold constant from last known value.
void LtiSdeGaussianProcessVectorSpacePrior::setExogenousInputs(const std::vector<Time>& times,
                                const std::vector<ValueType>& values){
  CHECK(false) << "Not implemented, todo Sean";
}

void LtiSdeGaussianProcessVectorSpacePrior::addKeyTime(const Time& time) {
  Coefficient coefficient;
  // integrate
  Key outKey = manager_.insertCoefficient(time, coefficient);
}

void LtiSdeGaussianProcessVectorSpacePrior::addKeyTimes(const std::vector<Time>& times) {
  std::vector<Key> outKeys;
  std::vector<Coefficient> coefficients(times.size());
  for (size_t i = 0; i < times.size(); ++i) {
    // integrate
  }
  manager_.insertCoefficients(times, coefficients, &outKeys);
}

Eigen::VectorXd LtiSdeGaussianProcessVectorSpacePrior::evaluate(Time time) const {
  return Eigen::VectorXd::Zero(this->dim(),1);
}

Eigen::VectorXd LtiSdeGaussianProcessVectorSpacePrior::evaluateDerivative(Time time, unsigned derivativeOrder) const {
  CHECK(false) << "Not implemented, todo Sean";
  return Eigen::VectorXd::Zero(this->dim(),1);
}

/// Evaluate the curve and interpolation matrices K(t)K^{-1}.
Eigen::VectorXd LtiSdeGaussianProcessVectorSpacePrior::evaluateAndInterpMatrices(Time time, const std::vector<Time>& keyTimes,
                                                                                 const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                                                 const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const {
  CHECK_EQ(keyTimes.size(), 2)                        << "Expected number of interpolation key times is not 2";
  CHECK_EQ(outEvalAtKeyTimes.size(), keyTimes.size()) << "Output vector is not initialized to have the correct number of entries";
  CHECK_EQ(outInterpMatrices.size(), keyTimes.size()) << "Output vector is not initialized to have the correct number of entries";

  CHECK_NOTNULL(outEvalAtKeyTimes[0]);
  CHECK_NOTNULL(outEvalAtKeyTimes[1]);
  *(outEvalAtKeyTimes[0]) = Eigen::VectorXd::Zero(this->dim(),1);
  *(outEvalAtKeyTimes[1]) = Eigen::VectorXd::Zero(this->dim(),1);

  CHECK_NOTNULL(outInterpMatrices[0]);
  CHECK_NOTNULL(outInterpMatrices[1]);
  *(outInterpMatrices[0]) = Eigen::MatrixXd::Identity(this->dim(),this->dim()) * 1.0;
  *(outInterpMatrices[1]) = Eigen::MatrixXd::Identity(this->dim(),this->dim()) * 0.0;

  return evaluate(time);
}

Eigen::VectorXd LtiSdeGaussianProcessVectorSpacePrior::evaluateDerivativeAndInterpMatrices(Time time, unsigned derivativeOrder, const std::vector<Time>& keyTimes,
                                                                                           const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                                                           const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const {
  CHECK(false) << "Not implemented, todo Sean";
  return Eigen::VectorXd::Zero(this->dim(),1);
}

/// \brief Get an evaluator at this time
LtiSdeGaussianProcessVectorSpacePrior::EvaluatorTypePtr LtiSdeGaussianProcessVectorSpacePrior::getEvaluator(const Time& time) const {
  CHECK(false) << "Gaussian Process Priors based on Stochastic Differential Equations do not have evaluators.";
  //boost::shared_ptr< Evaluator<VectorSpaceConfig> > rval( new Evaluator<VectorSpaceConfig>() );
  //return rval;
}

void LtiSdeGaussianProcessVectorSpacePrior::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

boost::unordered_map<Key, KeyCoefficientTime> LtiSdeGaussianProcessVectorSpacePrior::getKeyCoefficientTime() const {
  return manager_.getKeyCoefficientTime();
}

} // namespace curves
