#include <curves/LinearSdeGaussianProcessVectorSpacePrior.hpp>
#include <curves/LinearSdeGaussianProcessPriorFactorEvaluator.hpp>
#include <boost/make_shared.hpp>
#include <iostream>
#include <Eigen/Dense>

namespace curves {

LinearSdeGaussianProcessVectorSpacePrior::LinearSdeGaussianProcessVectorSpacePrior(size_t dimension, Eigen::MatrixXd stationaryPowerSpectralDensity) :
    GaussianProcessVectorSpacePrior(dimension), stationaryPowerSpectralDensity_(stationaryPowerSpectralDensity), invStationaryPowerSpectralDensity_(stationaryPowerSpectralDensity.inverse()) {
  CHECK_EQ(stationaryPowerSpectralDensity.rows(), stationaryPowerSpectralDensity.cols()) << "Stationary power spectral density matrix is not square.";
  initialized_ = false;
  initialTime_ = 0;
  initialMean_ = Eigen::VectorXd::Zero(dimension,1);
  initialCovariance_ = Eigen::MatrixXd::Zero(dimension,dimension);
}

LinearSdeGaussianProcessVectorSpacePrior::~LinearSdeGaussianProcessVectorSpacePrior() {}

void LinearSdeGaussianProcessVectorSpacePrior::print(const std::string& str) const {
  std::cout << "====================================" << std::endl;
  std::cout << "======= LINEAR SDE GP Prior ========" << std::endl;
  std::cout << str << std::endl;
  if (isInitialized()) {
    std::cout << "initial time: " << getMinTime() << std::endl;
    std::cout << "initial state: " << std::endl << initialMean_ << std::endl;
    std::cout << "initial cov: " << std::endl << initialCovariance_ << std::endl;
  } else {
    std::cout << "prior is uninitialized!" << std::endl;
  }
  std::cout << "dimension: " << dim() << std::endl;
  std::cout << "num of keytimes: " << getNumKeyTimes() << std::endl;
  std::cout << "power spectral density: " << std::endl << stationaryPowerSpectralDensity_ << std::endl;
  std::cout << "inv power spectral density: " << std::endl << invStationaryPowerSpectralDensity_ << std::endl;

}

Time LinearSdeGaussianProcessVectorSpacePrior::getMaxTime() const {
  CHECK(false) << "Not implemented, todo Sean";
  // return largest key time?
}

Time LinearSdeGaussianProcessVectorSpacePrior::getMinTime() const {
  return initialTime_;
}

void LinearSdeGaussianProcessVectorSpacePrior::initialize(Time initialTime, Eigen::VectorXd initialMean, Eigen::MatrixXd initialCovariance) {
  CHECK_EQ(initialMean.rows(), this->dim()) << "Dimension of initial mean does not match dimension of curve";
  CHECK_EQ(initialCovariance.rows(), this->dim()) << "Dimension of initial covariance does not match dimension of curve";
  CHECK_EQ(initialCovariance.rows(), initialCovariance.cols()) << "Initial covariance is not square.";
  initialized_ = true;
  initialTime_ = initialTime;
  initialMean_ = initialMean;
  initialCovariance_ = initialCovariance;
}

bool LinearSdeGaussianProcessVectorSpacePrior::isInitialized() const {
  return initialized_;
}

void LinearSdeGaussianProcessVectorSpacePrior::updateFromExogenousInputChange(Time time) {
  if (this->getNumKeyTimes() > 0)
  {
    // find where this input begins to affect keytime coefficients
    std::map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator it1 = keytimeToMeanSorted_.upper_bound(time);
    for(; it1 != keytimeToMeanSorted_.end(); ++it1)
    {
      // Get iterator to previous keyframe
      std::map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator it0 = it1; --it0;
      CHECK(it1->second->prevTime == it0->second->time) << "`Previous time' stored in keytime 1 does not match time of keytime 0.";

      // Update lifted exogenous inputs and mean function
      it1->second->liftedExogenousInput = calculateLiftedExogenousInput(it0->first, it1->first);
      it1->second->mean = it1->second->stateTransitionMatrix * it0->second->mean + it1->second->liftedExogenousInput;
    }
  }
}

void LinearSdeGaussianProcessVectorSpacePrior::addExogenousInput(Time time, const ValueType& value, bool updateMean) {

  // insert into timeToExogenousInputValue_
  timeToExogenousInputValue_.insert(std::pair<Time,Eigen::VectorXd>(time, value));

  // update if required
  if (updateMean) {
    this->updateFromExogenousInputChange(time);
  }
}

/// Append an exogenous input. Exogenous inputs create a step function that is assumed to hold constant from last known value.
void LinearSdeGaussianProcessVectorSpacePrior::addExogenousInput(Time time, const ValueType& value) {
  this->addExogenousInput(time, value, true);
}

/// Set the discrete time exogenous inputs.
void LinearSdeGaussianProcessVectorSpacePrior::setExogenousInputs(const std::vector<Time>& times, const std::vector<ValueType>& values){
  CHECK_EQ(times.size(), values.size()) << "times and values must be of the same size";
  CHECK_GT(times.size(), 0) << "vectors of times and values must have greater than length zero";
  timeToExogenousInputValue_.clear();
  for (size_t i = 0; i < values.size(); ++i) {
    this->addExogenousInput(times.at(i), values.at(i), false);
  }
  this->updateFromExogenousInputChange(keytimeToMeanSorted_.begin()->second->time);
}

void LinearSdeGaussianProcessVectorSpacePrior::addKeyTime(const Time& time, const Key& assocKey) {
  CHECK(isInitialized()) << "Prior has not yet been initialized.";
  CHECK_GE(time, initialTime_) << "New key time is less than or equal to the initial prior time.";

  if ( this->getNumKeyTimes() == 0 ) {
    if (time == initialTime_) {
      boost::shared_ptr<LinearSdeCoefficient> coeff = boost::shared_ptr<LinearSdeCoefficient> (new LinearSdeCoefficient());
      coeff->assocKey = assocKey;
      coeff->time = time;
      coeff->mean = initialMean_;
      coeff->inverseLiftedCovarianceMatrix = initialCovariance_.inverse();
      keytimeToMean_.insert(std::pair<Time,boost::shared_ptr<LinearSdeCoefficient> >(time, coeff));
      keytimeToMeanSorted_.insert(std::pair<Time,boost::shared_ptr<LinearSdeCoefficient> >(time, coeff));
    } else {
      CHECK(false) << "First keytime must match initialization time.";
      return;
    }
  } else {
    // Find insertion spot for time into data collection
    std::map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator it = keytimeToMeanSorted_.upper_bound(time);

    if(it == keytimeToMeanSorted_.end()) // insert new keytime at end
    {
      // Decrement iterator to last keytime entry
      --it;

      // Make new coefficient
      boost::shared_ptr<LinearSdeCoefficient> coeff = boost::shared_ptr<LinearSdeCoefficient> (new LinearSdeCoefficient());
      coeff->time = time;

      // Evaluate information relative to previous time
      coeff->prevTime = it->first;
      coeff->liftedExogenousInput = calculateLiftedExogenousInput(coeff->prevTime, coeff->time);
      coeff->stateTransitionMatrix = calculateStateTransitionMatrix(coeff->time, coeff->prevTime);
      coeff->inverseLiftedCovarianceMatrix = calculateInverseLiftedCovarianceMatrix(coeff->prevTime, coeff->time);

      // Evaluate mean function at new keytime
      coeff->mean = coeff->stateTransitionMatrix * it->second->mean + coeff->liftedExogenousInput;

      // Insert
      keytimeToMean_.insert(std::pair<Time,boost::shared_ptr<LinearSdeCoefficient> >(time, coeff));
      keytimeToMeanSorted_.insert(std::pair<Time,boost::shared_ptr<LinearSdeCoefficient> >(time, coeff));
    } else { // insert new keytime into middle
      CHECK(false) << "to be implemented.. lets try not to do this for now";
    }
  }
}

void LinearSdeGaussianProcessVectorSpacePrior::addKeyTimes(const std::vector<Time>& times, const std::vector<Key>& assocKeys) {
  CHECK_EQ(times.size(), assocKeys.size());
  for (size_t i = 0; i < times.size(); ++i) {
    this->addKeyTime(times.at(i), assocKeys.at(i));
  }
}

unsigned LinearSdeGaussianProcessVectorSpacePrior::getNumKeyTimes() const {
  return keytimeToMeanSorted_.size();
}

void LinearSdeGaussianProcessVectorSpacePrior::clearKeyTimes() {
  keytimeToMean_.clear();
  keytimeToMeanSorted_.clear();
}

std::vector<boost::shared_ptr<GaussianProcessPriorFactorEvaluator> > LinearSdeGaussianProcessVectorSpacePrior::getPriorFactors() const {

  std::vector<boost::shared_ptr<GaussianProcessPriorFactorEvaluator> > outFactors;
  std::map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator it = keytimeToMeanSorted_.begin();
  outFactors.push_back(boost::shared_ptr<GaussianProcessPriorFactorEvaluator>(new LinearSdeGaussianProcessPriorFactorEvaluator(it->second->assocKey, it->second)));
  Key lastKey = it->second->assocKey;
  for (; it != keytimeToMeanSorted_.end(); ++it) {
    outFactors.push_back(boost::shared_ptr<GaussianProcessPriorFactorEvaluator>(new LinearSdeGaussianProcessPriorFactorEvaluator(lastKey, it->second->assocKey, it->second)));
    lastKey = it->second->assocKey;
  }
  return outFactors;
}

Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluate(Time time) const {
  CHECK(isInitialized()) << "Prior has not yet been initialized.";
  CHECK_GE(time, initialTime_) << "New key time is less than or equal to the initial prior time.";
  CHECK_GT(getNumKeyTimes(), 0) << "Please add keytimes before using an SDE prior.";

  // See if we are evaluating at exactly a keytime
  boost::unordered_map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator keyIt = keytimeToMean_.find(time);
  if (keyIt != keytimeToMean_.end()) { // evaluate is at exactly a keytime
    return keyIt->second->mean;
  } else {
    // Get iterator to upper-bounding coefficient (in time)
    std::map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator it = keytimeToMeanSorted_.upper_bound(time);
    CHECK(it != keytimeToMeanSorted_.begin()) << "Tried to evaluate before initial time.";
    --it; // decrement iterator to lower-bounding coefficient
    return calculateStateTransitionMatrix(time, it->first) * it->second->mean + calculateLiftedExogenousInput(it->first, time);
  }

}

/// Evaluate the curve and interpolation matrices K(t)K^{-1}.
Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluateAndInterpMatrices(Time time, const std::vector<Time>& keyTimes,
                                                                                    const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                                                    const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const {
  CHECK_GT(getNumKeyTimes(), 0) << "Please add keytimes before using an SDE prior.";
  CHECK_EQ(keyTimes.size(), 2)  << "Expected number of interpolation key times is not 2";
  CHECK_GE(time, keyTimes[0])   << "Time is not in correct bounds";
  CHECK_LE(time, keyTimes[1])   << "Time is not in correct bounds";
  CHECK_EQ(outEvalAtKeyTimes.size(), keyTimes.size()) << "Output vector is not initialized to have the correct number of entries";
  CHECK_NOTNULL(outEvalAtKeyTimes[0]);
  CHECK_NOTNULL(outEvalAtKeyTimes[1]);
  CHECK_EQ(outInterpMatrices.size(), keyTimes.size()) << "Output vector is not initialized to have the correct number of entries";
  CHECK_NOTNULL(outInterpMatrices[0]);
  CHECK_NOTNULL(outInterpMatrices[1]);

  boost::unordered_map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator it0 = keytimeToMean_.find(keyTimes[0]);
  CHECK(it0 != keytimeToMean_.end()) << "Keytime did not exist in map";
  boost::unordered_map<Time, boost::shared_ptr<LinearSdeCoefficient> >::const_iterator it1 = keytimeToMean_.find(keyTimes[1]);
  CHECK(it1 != keytimeToMean_.end()) << "Keytime did not exist in map";
  CHECK(it1->second->prevTime == it0->second->time) << "`Previous time' stored in keytime 1 does not match time of keytime 0.";
  *(outEvalAtKeyTimes[0]) = it0->second->mean;
  *(outEvalAtKeyTimes[1]) = it1->second->mean;

  if (time == keyTimes[0]) {
    *(outInterpMatrices[0]) = Eigen::MatrixXd::Identity(this->dim(),this->dim());
    *(outInterpMatrices[1]) = Eigen::MatrixXd::Zero(this->dim(),this->dim());
  } else if (time == keyTimes[1]) {
    *(outInterpMatrices[0]) = Eigen::MatrixXd::Zero(this->dim(),this->dim());
    *(outInterpMatrices[1]) = Eigen::MatrixXd::Identity(this->dim(),this->dim());
  } else {
    *(outInterpMatrices[1]) = calculateLiftedCovarianceMatrix(it0->second->time, time) *
                              calculateStateTransitionMatrix(it1->second->time, time).transpose() *
                              it1->second->inverseLiftedCovarianceMatrix;
    *(outInterpMatrices[0]) = calculateStateTransitionMatrix(time, it0->second->time) - (*(outInterpMatrices[1])) * it1->second->stateTransitionMatrix;
  }

  return evaluate(time);
}

Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluateDerivative(Time time, unsigned derivativeOrder) const {
  CHECK_GT(getNumKeyTimes(), 0) << "Please add keytimes before using an SDE prior.";
  CHECK(false) << "Not implemented, todo Sean";
  return Eigen::VectorXd::Zero(this->dim(),1);
}

Eigen::VectorXd LinearSdeGaussianProcessVectorSpacePrior::evaluateDerivativeAndInterpMatrices(Time time, unsigned derivativeOrder, const std::vector<Time>& keyTimes,
                                                                                              const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                                                              const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const {
  CHECK_GT(getNumKeyTimes(), 0) << "Please add keytimes before using an SDE prior.";
  CHECK(false) << "Not implemented, todo Sean";
  return Eigen::VectorXd::Zero(this->dim(),1);
}

/// \brief Get an evaluator at this time
LinearSdeGaussianProcessVectorSpacePrior::EvaluatorTypePtr LinearSdeGaussianProcessVectorSpacePrior::getEvaluator(const Time& time) const {
  CHECK_GT(getNumKeyTimes(), 0) << "Please add keytimes before using an SDE prior.";
  CHECK(false) << "Gaussian Process Priors based on Stochastic Differential Equations do not have evaluators.";
  //boost::shared_ptr< Evaluator<VectorSpaceConfig> > rval( new Evaluator<VectorSpaceConfig>() );
  //return rval;
}

void LinearSdeGaussianProcessVectorSpacePrior::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

} // namespace curves
