/*
 * PolynomialSplineContainer.tpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */


namespace curves {


template <int splineOrder_>
PolynomialSplineContainer<splineOrder_>::PolynomialSplineContainer():
    timeOffset_(0.0),
    containerTime_(0.0),
    containerDuration_(0.0),
    activeSplineIdx_(0),
    hessian_(),
    linearTerm_(),
    equalityConstraintJacobian_(),
    equalityConstraintTargetValues_(),
    inequalityConstraintJacobian_(),
    inequalityConstraintMinValues_()
{
  // Make sure that the container is correctly emptied
  reset();

  minimizer_.reset(new numopt_quadprog::ActiveSetFunctionMinimizer());
  costFunction_.reset(new numopt_common::QuadraticObjectiveFunction());
  functionConstraints_.reset(new numopt_common::LinearFunctionConstraints());
  quadraticProblem_.reset(new numopt_common::QuadraticProblem(costFunction_, functionConstraints_));
}

template <int splineOrder_>
PolynomialSplineContainer<splineOrder_>::~PolynomialSplineContainer()
{

}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::advance(double dt)
{
  if (splines_.empty() || containerTime_ >= containerDuration_ || activeSplineIdx_ == splines_.size()) {
    return false;
  }

  // advance in time.
  containerTime_ += dt;

  // check if spline index needs to be increased.
  if ((containerTime_ - timeOffset_ >= splines_[activeSplineIdx_].getSplineDuration())) {
    if (activeSplineIdx_ < (splines_.size() - 1)) {
      timeOffset_ += splines_[activeSplineIdx_].getSplineDuration();
    }
    activeSplineIdx_++;
  }

  return true;
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::setContainerTime(double t)
{
  containerTime_ = t;
  double timeOffset;
  activeSplineIdx_ = getActiveSplineIndexAtTime(t, timeOffset);
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::reset()
{
  splines_.clear();
  activeSplineIdx_ = 0;
  containerDuration_ = 0.0;
  resetTime();
  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::resetTime()
{
  timeOffset_ = 0.0;
  containerTime_ = 0.0;
  activeSplineIdx_ = 0;
  return true;
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getContainerDuration() const
{
  return containerDuration_;
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getContainerTime() const
{
  return containerTime_;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::isEmpty() const
{
  return splines_.empty();
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getPosition() const
{
  if (splines_.empty()) {
    return 0.0;
  }
  if (activeSplineIdx_ == splines_.size()) {
    return splines_.back().getPositionAtTime(containerTime_ - timeOffset_);
  }
  return splines_[activeSplineIdx_].getPositionAtTime(containerTime_ - timeOffset_);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getVelocity() const {
  if (splines_.empty()) {
    return 0.0;
  }
  if (activeSplineIdx_ == splines_.size()) {
    return splines_.back().getVelocityAtTime(containerTime_ - timeOffset_);
  }
  return splines_[activeSplineIdx_].getVelocityAtTime(containerTime_ - timeOffset_);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getAcceleration() const
{
  if (splines_.empty()) {
    return 0.0;
  }
  if (activeSplineIdx_ == splines_.size()) {
    return splines_.back().getAccelerationAtTime(containerTime_ - timeOffset_);
  }
  return splines_[activeSplineIdx_].getAccelerationAtTime(containerTime_ - timeOffset_);
}

template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getActiveSplineIndexAtTime(double t, double& timeOffset) const {
  if (splines_.empty()) return -1;
  timeOffset = 0.0;

  for (size_t i = 0; i < splines_.size(); i++) {
    if ((t - timeOffset < splines_[i].getSplineDuration())) {
      return i;
    }
    if (i < (splines_.size() - 1)) {
      timeOffset += splines_[i].getSplineDuration();
    }
  }

  return (splines_.size() - 1);
}

template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getActiveSplineIndex() const {
  return activeSplineIdx_;
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getPositionAtTime(double t) const {
  double timeOffset = 0.0;
  const int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);

  if (activeSplineIdx < 0) {
    // Spline container is empty.
    return 0.0;
  }

  if (activeSplineIdx == splines_.size()) {
    return splines_.back().getPositionAtTime(splines_.back().getSplineDuration());
  }

  return splines_[activeSplineIdx].getPositionAtTime(t - timeOffset);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getVelocityAtTime(double t) const
{
  double timeOffset = 0.0;
  const int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);

  if (activeSplineIdx < 0) {
    // Spline container is empty.
    return 0.0;
  }

  if (activeSplineIdx == splines_.size()) {
    return splines_.back().getVelocityAtTime(splines_.back().getSplineDuration());
  }

  return splines_[activeSplineIdx].getVelocityAtTime(t - timeOffset);
}


template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getAccelerationAtTime(double t) const
{
  double timeOffset = 0.0;
  const int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);

  if (activeSplineIdx < 0) {
    // Spline container is empty.
    return 0.0;
  }

  if (activeSplineIdx == splines_.size()) {
    return splines_.back().getAccelerationAtTime(splines_.back().getSplineDuration());
  }

  return splines_[activeSplineIdx].getAccelerationAtTime(t - timeOffset);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getEndPosition() const {
  if (splines_.empty()) {
    // Spline container is empty.
    return 0.0;
  }

  const double lastSplineDuration = splines_.back().getSplineDuration();
  return splines_.back().getPositionAtTime(lastSplineDuration);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getEndVelocity() const {
  if (splines_.empty()) {
    // Spline container is empty.
    return 0.0;
  }

  const double lastSplineDuration = splines_.back().getSplineDuration();
  return splines_.back().getVelocityAtTime(lastSplineDuration);
}

template <int splineOrder_>
double PolynomialSplineContainer<splineOrder_>::getEndAcceleration() const {
  if (splines_.empty()) {
    // Spline container is empty.
    return 0.0;
  }

  const double lastSplineDuration = splines_.back().getSplineDuration();
  return splines_.back().getAccelerationAtTime(lastSplineDuration);
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setData(const std::vector<double>& knotDurations,
             const std::vector<double>& knotPositions,
             double initialVelocity,
             double initialAcceleration,
             double finalVelocity,
             double finalAcceleration,
             double weightMinAccel) {

  bool success = true;

  success &= reset();

  // Set up optimization parameters
  const unsigned int numSplines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  const unsigned int solutionSpaceDimension = numSplines*num_coeffs_spline;
  const unsigned int num_junctions = numSplines-1;

  if (numSplines<1) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Not sufficient knot points are available!");
    return false;
  }

  // total number of constraints
  constexpr unsigned int num_initial_constraints = 3;   // pos, vel, accel
  constexpr unsigned int num_final_constraints = 3;     // pos, vel, accel
  constexpr unsigned int num_constraint_junction = 4;   // pos (2x), vel, accel
  const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
  const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

  // drop constraints if necessary
  if (num_constraints>solutionSpaceDimension) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Number of equality constraints is larger than number of coefficients. Drop acceleration constraints!");
    return setData(knotDurations, knotPositions, initialVelocity, finalVelocity, weightMinAccel);
  }

  // Vector containing durations of splines
  std::vector<double> splineDurations(numSplines);
  for (unsigned int splineId=0; splineId<numSplines; splineId++) {
    splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

    if (splineDurations[splineId]<=0.0) {
      MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Invalid spline duration at index" << splineId << ": " << splineDurations[splineId]);
      return false;
    }
  }

  // Initialize Equality matrices
  equalityConstraintJacobian_.setZero(num_constraints, solutionSpaceDimension);
  equalityConstraintTargetValues_.setZero(num_constraints);
  unsigned int constraintIdx = 0;


  // Initial conditions
  addInitialConditions(
      (Eigen::VectorXd() << knotPositions.front(), initialVelocity, initialAcceleration).finished(),
      constraintIdx);

  // Final conditions
  addFinalConditions(
      (Eigen::VectorXd() << knotPositions.back(), finalVelocity, finalAcceleration).finished(),
      constraintIdx,  splineDurations.back(), num_junctions);

  // Junction conditions
  addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

  if (num_constraints != constraintIdx) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Wrong number of equality constraints!");
    return false;
  }

  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(solutionSpaceDimension);

  // Minimize spline coefficients
  if (weightMinAccel<=0.0 || num_constraints==solutionSpaceDimension) {
    coeffs = equalityConstraintJacobian_.colPivHouseholderQr().solve(equalityConstraintTargetValues_);
  }

  // Minimize acceleration
  else {
    success &= addObjective(splineDurations, solutionSpaceDimension, numSplines, weightMinAccel);
    success &= setUpOptimizationMatrices(solutionSpaceDimension);
    success &= solveQP(coeffs);
  }

  // Extract spline coefficients and add splines
  success &= extractSplineCoefficients(coeffs, splineDurations, numSplines);

  return success;

}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setData(const std::vector<double>& knotDurations,
             const std::vector<double>& knotPositions,
             double initialVelocity,
             double finalVelocity,
             double weightMinAccel) {

  bool success = true;

  success &= reset();

  // Set up optimization parameters
  const unsigned int num_splines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  const unsigned int num_coeffs = num_splines*num_coeffs_spline;
  const unsigned int num_junctions = num_splines-1;

  if (num_splines<1) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Not sufficient knot points are available!");
    return false;
  }

  // total number of constraints
  constexpr unsigned int num_initial_constraints = 2;   // pos, vel
  constexpr unsigned int num_final_constraints = 2;     // pos, vel
  constexpr unsigned int num_constraint_junction = 3;   // pos (2x), vel
  const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
  const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

  // drop constraints if necessary
  if (num_constraints>num_coeffs) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Number of equality constraints is larger than number of coefficients. Drop velocity constraints!");
    return setData(knotDurations, knotPositions);
  }

  // Vector containing durations of splines
  std::vector<double> splineDurations(num_splines);
  for (unsigned int splineId=0; splineId<num_splines; splineId++) {
    splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

    if (splineDurations[splineId]<=0.0) {
      MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Invalid spline duration at index" << splineId << ": " << splineDurations[splineId]);
      return false;
    }
  }

  // Initialize Equality matrices
  equalityConstraintJacobian_.setZero(num_constraints, num_coeffs);
  equalityConstraintTargetValues_.setZero(num_constraints);
  unsigned int constraintIdx = 0;

  // Initial conditions
  addInitialConditions(
      (Eigen::VectorXd() << knotPositions.front(), initialVelocity).finished(),
      constraintIdx);

  // Final conditions
  addFinalConditions(
      (Eigen::VectorXd() << knotPositions.back(), finalVelocity).finished(),
      constraintIdx,  splineDurations.back(), num_junctions);

  // Junction conditions
  addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

  if (num_constraints!=constraintIdx) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData] Wrong number of equality constraints!");
    return false;
  }

  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(num_coeffs);

  // Minimize spline coefficients
  if (weightMinAccel<=0.0 || num_constraints==num_coeffs) {
    coeffs = equalityConstraintJacobian_.colPivHouseholderQr().solve(equalityConstraintTargetValues_);
  }

  // Minimize acceleration
  else {
    success &= addObjective(splineDurations, num_coeffs, num_splines, weightMinAccel);
    success &= setUpOptimizationMatrices(num_coeffs);
    success &= solveQP(coeffs);
  }

  // Extract spline coefficients and add splines
  success &= extractSplineCoefficients(coeffs, splineDurations, num_splines);

  return success;

}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setData(const std::vector<double>& knotDurations,
             const std::vector<double>& knotPositions) {

  bool success = true;
  const int num_splines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  typename SplineType::SplineCoefficients coefficients;

  if (splineOrder_==0 || num_splines<1) {
    return false;
  }

  std::fill(coefficients.begin(), coefficients.end(), 0.0);

  success &= reset();

  for (unsigned int splineId = 0; splineId<num_splines; splineId++) {
    const double duration = knotDurations[splineId+1]-knotDurations[splineId];

    if (duration<=0.0) {
      return false;
    }

    coefficients[num_coeffs_spline-1] = knotPositions[splineId]; //a0
    coefficients[num_coeffs_spline-2] = (knotPositions[splineId+1]-knotPositions[splineId]) / duration; // a1
    success &= this->addSpline(SplineType(coefficients, duration));
  }

  return success;

}


template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getCoeffIndex(int splineIdx, int aIdx) const {
  return splineIdx*(splineOrder_+1) + aIdx;
}

template <int splineOrder_>
int PolynomialSplineContainer<splineOrder_>::getSplineColumnIndex(int splineIdx) const
{
  return getCoeffIndex(splineIdx, 0);
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::addInitialConditions(const Eigen::VectorXd& initialConditions,
                          unsigned int& constraintIdx) {

  // time container
  typename SplineType::EigenTimeVectorType timeVec;

  // initial position
  if (initialConditions.size()>0) {
    SplineType::getTimeVector(timeVec, 0.0);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = initialConditions(0);
    constraintIdx++;
  }

  // initial velocity
  if (initialConditions.size()>1) {
    SplineType::getdTimeVector(timeVec, 0.0);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = initialConditions(1);
    constraintIdx++;
  }

  // initial acceleration
  if (initialConditions.size()>2) {
    SplineType::getddTimeVector(timeVec, 0.0);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = initialConditions(2);
    constraintIdx++;
  }
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::addFinalConditions(const Eigen::VectorXd& finalConditions,
                        unsigned int& constraintIdx,
                        double lastSplineDuration,
                        unsigned int lastSplineId) {
  // time container
  typename SplineType::EigenTimeVectorType timeVec;

  // initial position
  if (finalConditions.size()>0) {
    SplineType::getTimeVector(timeVec, lastSplineDuration);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = finalConditions(0);
    constraintIdx++;
  }

  // initial velocity
  if (finalConditions.size()>1) {
    SplineType::getdTimeVector(timeVec, lastSplineDuration);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = finalConditions(1);
    constraintIdx++;
  }

  // initial acceleration
  if (finalConditions.size()>2) {
    SplineType::getddTimeVector(timeVec, lastSplineDuration);
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec;
    equalityConstraintTargetValues_(constraintIdx) = finalConditions(2);
    constraintIdx++;
  }
}

template <int splineOrder_>
void PolynomialSplineContainer<splineOrder_>::addJunctionsConditions(const std::vector<double>& splineDurations,
                            const std::vector<double>& knotPositions,
                            unsigned int& constraintIdx,
                            unsigned int num_junctions) {

  // time containers
  typename SplineType::EigenTimeVectorType timeVec0, dTimeVec0, ddTimeVec0;
  typename SplineType::EigenTimeVectorType timeVecTf, dTimeVecTf, ddTimeVecTf;

  // get time container at zero time
  SplineType::getTimeVectorAtZero(timeVec0);
  SplineType::getdTimeVectorAtZero(dTimeVec0);
  SplineType::getddTimeVectorAtZero(ddTimeVec0);


  for (unsigned int splineId=0; splineId<num_junctions; splineId++) {
    const unsigned int nextSplineId = splineId+1;

    // get time container at spline duration
    SplineType::getTimeVector(timeVecTf, splineDurations[splineId]);
    SplineType::getdTimeVector(dTimeVecTf, splineDurations[splineId]);
    SplineType::getddTimeVector(ddTimeVecTf, splineDurations[splineId]);

    // smooth position transition with fixed positions
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  timeVecTf;
    equalityConstraintTargetValues_(constraintIdx) = knotPositions[nextSplineId];
    constraintIdx++;

    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) =  timeVec0;
    equalityConstraintTargetValues_(constraintIdx) = knotPositions[nextSplineId];
    constraintIdx++;

    // smooth velocity transition
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  dTimeVecTf;
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) = -dTimeVec0;
    equalityConstraintTargetValues_(constraintIdx) = 0.0;
    constraintIdx++;

    // smooth acceleration transition
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  ddTimeVecTf;
    equalityConstraintJacobian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) = -ddTimeVec0;
    equalityConstraintTargetValues_(constraintIdx) = 0.0;
    constraintIdx++;
  }
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::getAccelerationMinimizerBlock(MinAccMat& mat, double tf) const {
  if (splineOrder_ == 5) {
    const double tf2 = boost::math::pow<2>(tf);
    const double tf3 = tf2*tf;
    const double tf4 = tf3*tf;
    const double tf5 = tf4*tf;
    const double tf6 = tf5*tf;
    const double tf7 = tf6*tf;

    mat << 400.0/7.0*tf7, 40.0*tf6,       24.0*tf5,       10.0*tf4,
           40.0*tf6,      28.8*tf5,       18.0*tf4,       8.0*tf3,
           24.0*tf5,      18.0*tf4,       12.0*tf3,       6.0*tf2,
           10.0*tf4,      8.0*tf3,        6.0*tf2,        4.0*tf;
  }

  else if (splineOrder_ == 4) {
    const double tf2 = boost::math::pow<2>(tf);
    const double tf3 = tf2*tf;
    const double tf4 = tf3*tf;
    const double tf5 = tf4*tf;

    mat << 28.8*tf5,       18.0*tf4,       8.0*tf3,
           18.0*tf4,       12.0*tf3,       6.0*tf2,
           8.0*tf3,        6.0*tf2,        4.0*tf;

  }

  else if (splineOrder_ == 3) {
    const double tf2 = boost::math::pow<2>(tf);
    const double tf3 = tf2*tf;

    mat << 12.0*tf3,       6.0*tf2,
           6.0*tf2,        4.0*tf;
  }

  else if (splineOrder_ == 2) {
    mat << 4.0*tf;
  }

  else {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData::getAccelerationMinimizerBlock] Function has not been implemented so far.");
    return false;
  }

  return true;
}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::addObjective(
    const std::vector<double>& splineDurations,
    unsigned int num_coeffs,
    unsigned int num_splines,
    double weightMinAccel) {

  bool success = true;

  // Objective -> minimize acceleration along trajectory
  hessian_ = Eigen::MatrixXd::Zero(num_coeffs,num_coeffs);
  linearTerm_ = Eigen::VectorXd::Zero(num_coeffs);

  for (unsigned int splineId=0; splineId<num_splines; splineId++) {

    // Minimize acceleration
    MinAccMat coreMatrix;
    success &= getAccelerationMinimizerBlock(coreMatrix, splineDurations[splineId]);
    hessian_.block<SplineType::coefficientCount-2, SplineType::coefficientCount-2>(
        getSplineColumnIndex(splineId),  getSplineColumnIndex(splineId)) =
            coreMatrix*weightMinAccel;

    // Regularization for position and velocity -> cfMat must be positive definite.
    const unsigned int accelIdx = getCoeffIndex(splineId, SplineType::coefficientCount-2);
    hessian_.block<2,2>(accelIdx, accelIdx) += 1e-7*Eigen::Matrix2d::Identity();
  }

  return success;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::setUpOptimizationMatrices(unsigned int num_coeffs) {
  // Inequality constraints -> We don't use them
  inequalityConstraintJacobian_.setZero(0, num_coeffs);
  inequalityConstraintMinValues_.setZero(0);

  // Add quadratic problem
  costFunction_->setGlobalHessian(hessian_.sparseView());
  costFunction_->setLinearTerm(linearTerm_);
  functionConstraints_->setGlobalEqualityConstraintJacobian(equalityConstraintJacobian_.sparseView());
  functionConstraints_->setEqualityConstraintTargetValues(equalityConstraintTargetValues_);
  functionConstraints_->setGlobalInequalityConstraintJacobian(inequalityConstraintJacobian_.sparseView());
  functionConstraints_->setInequalityConstraintMinValues(inequalityConstraintMinValues_);
  functionConstraints_->setInequalityConstraintMaxValues(inequalityConstraintMinValues_);

  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::solveQP(Eigen::VectorXd& coeffs) {
  bool success = true;

  double cost = 0.0;
  numopt_common::ParameterizationIdentity params(coeffs);
  success &= minimizer_->minimize(quadraticProblem_.get(), params, cost);

  if (!success) {
    MELO_WARN_STREAM("[PolynomialSplineContainer::setData::solveQP] Failed to solve optimization!");
  } else {
    coeffs = params.getParams();
  }

  return success;
}


template <int splineOrder_>
bool PolynomialSplineContainer<splineOrder_>::extractSplineCoefficients(
    const Eigen::VectorXd& coeffs,
    const std::vector<double>& splineDurations,
    const unsigned int numSplines) {
  typename SplineType::SplineCoefficients coefficients;

  for (unsigned int splineId = 0; splineId<numSplines; ++splineId) {
    Eigen::Map<Eigen::VectorXd>(coefficients.data(), SplineType::coefficientCount, 1) =
        coeffs.segment<SplineType::coefficientCount>(getSplineColumnIndex(splineId));
    this->addSpline(SplineType(coefficients,splineDurations[splineId]));
  }

  return true;
}


} /* namespace */
