#include <curves/LinearSdeGaussianProcessVectorSpacePrior.hpp>
#include <boost/make_shared.hpp>
#include <iostream>

namespace curves {

LinearSdeGaussianProcessVectorSpacePrior::LinearSdeGaussianProcessVectorSpacePrior(size_t dimension, Eigen::MatrixXd stationaryPowerSpectralDensity) :
    GaussianProcessVectorSpacePrior(dimension), stationaryPowerSpectralDensity_(stationaryPowerSpectralDensity) {
  CHECK_EQ(stationaryPowerSpectralDensity.rows(), stationaryPowerSpectralDensity.cols()) << "Stationary power spectral density matrix is not square.";
}

LinearSdeGaussianProcessVectorSpacePrior::~LinearSdeGaussianProcessVectorSpacePrior() {}

void LinearSdeGaussianProcessVectorSpacePrior::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=======GAUSSIAN PROCESS PRIOR========" << std::endl;
  std::cout << "=========================================" <<std::endl;
}

Time LinearSdeGaussianProcessVectorSpacePrior::getMaxTime() const {
  CHECK(false) << "Not implemented, todo Sean";
  // return largest key time?
}

Time LinearSdeGaussianProcessVectorSpacePrior::getMinTime() const {
  return initialTime_;
}

void LinearSdeGaussianProcessVectorSpacePrior::initialize(Time initialTime, Eigen::VectorXd initialMean,
                                                          Eigen::MatrixXd initialInverseCovariance) {
  initialTime_ = initialTime;
  initialInverseCovariance_ = initialInverseCovariance;
  boost::shared_ptr<Eigen::VectorXd> sharedInitialMean(new Eigen::VectorXd(initialMean));
  keytimeToMeanEvaluation_.insert(keytimeToMeanEvaluation_.begin(), std::pair<Time, boost::shared_ptr<Eigen::VectorXd> >(initialTime, sharedInitialMean));
}

/// Append an exogenous input.
/// Exogenous inputs create a step function that is assumed to hold constant from last known value.
void LinearSdeGaussianProcessVectorSpacePrior::addExogenousInput(Time time, const ValueType& value) {
  CHECK(false) << "Not implemented, todo Sean";
}

/// Set the discrete time exogenous inputs.
/// Exogenous inputs create a step function that is assumed to hold constant from last known value.
void LinearSdeGaussianProcessVectorSpacePrior::setExogenousInputs(const std::vector<Time>& times,
                                const std::vector<ValueType>& values){
  CHECK(false) << "Not implemented, todo Sean";
}

void LinearSdeGaussianProcessVectorSpacePrior::addKeyTime(const Time& time) {
  CHECK(!keytimeToMeanEvaluation_.empty()) << "Prior has not yet been initialized.";
  CHECK_GT(time, initialTime_) << "New key time is less than or equal to the initial prior time.";

  // find insertion spot for time into keytimeToMeanEvaluation_
  std::map<Time, boost::shared_ptr<Eigen::VectorXd> >::const_iterator it = keytimeToMeanEvaluation_.upper_bound(time);
  CHECK(it != keytimeToMeanEvaluation_.begin()) << "Tried to insert keytime before initial time.";

  if(it == keytimeToMeanEvaluation_.end()) // insert new keytime at end
  {
    --it; // decrement iterator to last keytime entry

    // Evaluate state transition matrix and integrated exogenous inputs
    boost::shared_ptr<Eigen::MatrixXd> stateTransition = boost::make_shared<Eigen::MatrixXd> (calculateStateTransitionMatrix(time, it->first));
    boost::shared_ptr<Eigen::VectorXd> liftedExoInput = boost::make_shared<Eigen::VectorXd> (calculateLiftedExogenousInput(it->first, time));

    // Evaluate mean function at new keytime
    Eigen::VectorXd asdf = (*stateTransition) * (*it->second) + (*liftedExoInput);
    boost::shared_ptr<Eigen::VectorXd> meanEval = boost::make_shared<Eigen::VectorXd> ( asdf );

    // Evaluate inverse of lifted covariance matrix
    boost::shared_ptr<Eigen::MatrixXd> invLiftedCov = boost::make_shared<Eigen::MatrixXd> (calculateInverseLiftedCovarianceMatrix(it->first, time));

    // Create pair key
    std::pair<Time,Time> keytimePair(it->first, time);

    // Insert into maps
    keytimeToMeanEvaluation_.insert(keytimeToMeanEvaluation_.end(), std::pair<Time,boost::shared_ptr<Eigen::VectorXd> >(time, meanEval));
    keytimesToLiftedExogenousInput_.insert(std::pair<std::pair<Time,Time>, boost::shared_ptr<Eigen::VectorXd> >(keytimePair, liftedExoInput));
    keytimesToStateTransitionMatrix_.insert(std::pair<std::pair<Time,Time>, boost::shared_ptr<Eigen::MatrixXd> >(keytimePair, stateTransition));
    keytimesToInverseLiftedCovarianceMatrix_.insert(std::pair<std::pair<Time,Time>, boost::shared_ptr<Eigen::MatrixXd> >(keytimePair, invLiftedCov));

  } else { // insert new keytime into middle

    CHECK(false) << "to be implemented.. lets try not to do this for now";

  }
}

void LinearSdeGaussianProcessVectorSpacePrior::addKeyTimes(const std::vector<Time>& times) {
  for (size_t i = 0; i < times.size(); ++i) {
    this->addKeyTime(times.at(i));
  }
}

Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluate(Time time) const {
  std::map<Time, boost::shared_ptr<Eigen::VectorXd> >::const_iterator it = keytimeToMeanEvaluation_.upper_bound(time);
  CHECK(it != keytimeToMeanEvaluation_.begin()) << "Tried to insert keytime before initial time.";
  --it;
  return calculateStateTransitionMatrix(time, it->first)*(*it->second) + calculateLiftedExogenousInput(it->first, time);
}

Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluateDerivative(Time time, unsigned derivativeOrder) const {
  CHECK(false) << "Not implemented, todo Sean";
  return Eigen::VectorXd::Zero(this->dim(),1);
}

/// Evaluate the curve and interpolation matrices K(t)K^{-1}.
Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluateAndInterpMatrices(Time time, const std::vector<Time>& keyTimes,
                                                                                    const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                                                    const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const {
  CHECK_EQ(keyTimes.size(), 2)                        << "Expected number of interpolation key times is not 2";
  CHECK_EQ(outEvalAtKeyTimes.size(), keyTimes.size()) << "Output vector is not initialized to have the correct number of entries";
  CHECK_EQ(outInterpMatrices.size(), keyTimes.size()) << "Output vector is not initialized to have the correct number of entries";

  /// \todo should check that keyTimes exist... but log(n) search is a little slow...

  CHECK_NOTNULL(outEvalAtKeyTimes[0]);
  CHECK_NOTNULL(outEvalAtKeyTimes[1]);
  *(outEvalAtKeyTimes[0]) = *keytimeToMeanEvaluation_.at(keyTimes[0]);
  *(outEvalAtKeyTimes[1]) = *keytimeToMeanEvaluation_.at(keyTimes[1]);

  CHECK_NOTNULL(outInterpMatrices[0]);
  CHECK_NOTNULL(outInterpMatrices[1]);
  std::pair<Time,Time> keytimePair(keyTimes[0],keyTimes[1]);
  *(outInterpMatrices[1]) = calculateLiftedCovarianceMatrix(keyTimes[0],time) *
                            calculateStateTransitionMatrix(keyTimes[1], time).transpose() *
                            (*keytimesToInverseLiftedCovarianceMatrix_.at(keytimePair));
  *(outInterpMatrices[0]) = calculateStateTransitionMatrix(time, keyTimes[0]) - (*(outInterpMatrices[1])) * (*keytimesToStateTransitionMatrix_.at(keytimePair));

  return evaluate(time);
}

Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluateDerivativeAndInterpMatrices(Time time, unsigned derivativeOrder, const std::vector<Time>& keyTimes,
                                                                                              const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                                                              const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const {
  CHECK(false) << "Not implemented, todo Sean";
  return Eigen::VectorXd::Zero(this->dim(),1);
}

/// \brief Get an evaluator at this time
LinearSdeGaussianProcessVectorSpacePrior::EvaluatorTypePtr LinearSdeGaussianProcessVectorSpacePrior::getEvaluator(const Time& time) const {
  CHECK(false) << "Gaussian Process Priors based on Stochastic Differential Equations do not have evaluators.";
  //boost::shared_ptr< Evaluator<VectorSpaceConfig> > rval( new Evaluator<VectorSpaceConfig>() );
  //return rval;
}

void LinearSdeGaussianProcessVectorSpacePrior::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

} // namespace curves
