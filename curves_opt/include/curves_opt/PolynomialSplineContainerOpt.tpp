/*
 * PolynomialSplineContainerOpt.tpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */


namespace curves {


template <int splineOrder_>
PolynomialSplineContainerOpt<splineOrder_>::PolynomialSplineContainerOpt():
    Base(),
    hessian_(),
    linearTerm_(),
    inequalityConstraintJacobian_(),
    inequalityConstraintMinValues_()
{
  // Make sure that the container is correctly emptied.
  this->reset();

  minimizer_.reset(new numopt_quadprog::ActiveSetFunctionMinimizer());
  costFunction_.reset(new numopt_common::QuadraticObjectiveFunction());
  functionConstraints_.reset(new numopt_common::LinearFunctionConstraints());
  quadraticProblem_.reset(new numopt_common::QuadraticProblem(costFunction_, functionConstraints_));
}

template <int splineOrder_>
PolynomialSplineContainerOpt<splineOrder_>::~PolynomialSplineContainerOpt()
{

}

template <int splineOrder_>
bool PolynomialSplineContainerOpt<splineOrder_>::setDataOptimized(
    const std::vector<double>& knotDurations,
    const std::vector<double>& knotPositions,
    double initialVelocity, double initialAcceleration,
    double finalVelocity, double finalAcceleration,
    double weightMinAccel) {

  bool success = true;

  success &= this->reset();

  // Set up optimization parameters.
  const unsigned int numSplines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  const unsigned int solutionSpaceDimension = numSplines*num_coeffs_spline;
  const unsigned int num_junctions = numSplines-1;

  if (numSplines<1) {
    std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Not sufficient knot points are available!" << std::endl;
    return false;
  }

  // total number of constraints.
  constexpr unsigned int num_initial_constraints = 3;   // pos, vel, accel
  constexpr unsigned int num_final_constraints = 3;     // pos, vel, accel
  constexpr unsigned int num_constraint_junction = 4;   // pos (2x), vel, accel
  const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
  const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

  // Drop constraints if necessary.
  if (num_constraints>solutionSpaceDimension) {
    std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Number of equality constraints is larger than number of coefficients. Drop acceleration constraints!" << std::endl;
    return setDataOptimized(knotDurations, knotPositions, initialVelocity, finalVelocity, weightMinAccel);
  }

  // Vector containing durations of splines.
  std::vector<double> splineDurations(numSplines);
  for (unsigned int splineId=0; splineId<numSplines; splineId++) {
    splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

    if (splineDurations[splineId]<=0.0) {
      std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Invalid spline duration at index" << splineId << ": " << splineDurations[splineId] << std::endl;
      return false;
    }
  }

  // Initialize Equality matrices.
  this->equalityConstraintJacobian_.setZero(num_constraints, solutionSpaceDimension);
  this->equalityConstraintTargetValues_.setZero(num_constraints);
  unsigned int constraintIdx = 0;


  // Initial conditions.
  Eigen::VectorXd initialConditions(3);
  initialConditions << knotPositions.front(), initialVelocity, initialAcceleration;
  this->addInitialConditions(initialConditions, constraintIdx);

  // Final conditions.
  Eigen::VectorXd finalConditions(3);
  finalConditions << knotPositions.back(), finalVelocity, finalAcceleration;
  this->addFinalConditions(finalConditions, constraintIdx, splineDurations.back(), num_junctions);

  // Junction conditions.
  this->addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

  if (num_constraints != constraintIdx) {
    std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Wrong number of equality constraints!" << std::endl;
    return false;
  }

  Eigen::VectorXd coeffs;

  if (weightMinAccel<=0.0 || num_constraints==solutionSpaceDimension) {
    // Minimize spline coefficients.
    coeffs = this->equalityConstraintJacobian_.colPivHouseholderQr().solve(this->equalityConstraintTargetValues_);
  } else {
    // Minimize acceleration.
    success &= addObjective(splineDurations, solutionSpaceDimension, numSplines, weightMinAccel);
    success &= setUpOptimizationMatrices(solutionSpaceDimension);
    success &= solveQP(coeffs);
  }

  // Extract spline coefficients and add splines.
  success &= this->extractSplineCoefficients(coeffs, splineDurations, numSplines);

  return success;

}


template <int splineOrder_>
bool PolynomialSplineContainerOpt<splineOrder_>::setDataOptimized(
    const std::vector<double>& knotDurations,
    const std::vector<double>& knotPositions,
    double initialVelocity, double finalVelocity,
    double weightMinAccel) {

  bool success = true;

  success &= this->reset();

  // Set up optimization parameters.
  const unsigned int num_splines = knotDurations.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  const unsigned int num_coeffs = num_splines*num_coeffs_spline;
  const unsigned int num_junctions = num_splines-1;

  if (num_splines<1) {
    std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Not sufficient knot points are available!" << std::endl;
    return false;
  }

  // Total number of constraints.
  constexpr unsigned int num_initial_constraints = 2;   // pos, vel
  constexpr unsigned int num_final_constraints = 2;     // pos, vel
  constexpr unsigned int num_constraint_junction = 3;   // pos (2x), vel
  const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
  const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

  // Drop constraints if necessary.
  if (num_constraints>num_coeffs) {
    std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Number of equality constraints is larger than number of coefficients. Drop velocity constraints!" << std::endl;
    return this->setData(knotDurations, knotPositions);
  }

  // Vector containing durations of splines.
  std::vector<double> splineDurations(num_splines);
  for (unsigned int splineId=0; splineId<num_splines; splineId++) {
    splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

    if (splineDurations[splineId]<=0.0) {
      std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Invalid spline duration at index" << splineId << ": " << splineDurations[splineId] << std::endl;
      return false;
    }
  }

  // Initialize Equality matrices.
  this->equalityConstraintJacobian_.setZero(num_constraints, num_coeffs);
  this->equalityConstraintTargetValues_.setZero(num_constraints);
  unsigned int constraintIdx = 0;

  // Initial conditions.
  this->addInitialConditions(
      (Eigen::VectorXd() << knotPositions.front(), initialVelocity).finished(),
      constraintIdx);

  // Final conditions.
  this->addFinalConditions(
      (Eigen::VectorXd() << knotPositions.back(), finalVelocity).finished(),
      constraintIdx,  splineDurations.back(), num_junctions);

  // Junction conditions.
  this->addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

  if (num_constraints!=constraintIdx) {
    std::cout << "[PolynomialSplineContainerOpt::setDataOptimized] Wrong number of equality constraints!" << std::endl;
    return false;
  }

  Eigen::VectorXd coeffs;

  if (weightMinAccel<=0.0 || num_constraints==num_coeffs) {
    // Minimize spline coefficients.
    coeffs = this->equalityConstraintJacobian_.colPivHouseholderQr().solve(this->equalityConstraintTargetValues_);
  } else {
    // Minimize acceleration.
    success &= addObjective(splineDurations, num_coeffs, num_splines, weightMinAccel);
    success &= setUpOptimizationMatrices(num_coeffs);
    success &= solveQP(coeffs);
  }

  // Extract spline coefficients and add splines.
  success &= this->extractSplineCoefficients(coeffs, splineDurations, num_splines);

  return success;

}

template <int splineOrder_>
bool PolynomialSplineContainerOpt<splineOrder_>::addObjective(
    const std::vector<double>& splineDurations,
    unsigned int num_coeffs,
    unsigned int num_splines,
    double weightMinAccel) {

  bool success = true;

  // Minimize acceleration along trajectory.
  hessian_ = Eigen::MatrixXd::Zero(num_coeffs,num_coeffs);
  linearTerm_ = Eigen::VectorXd::Zero(num_coeffs);

  for (unsigned int splineId=0; splineId<num_splines; ++splineId) {
    // Minimize acceleration.
    MinAccMat coreMatrix;
    success &= getAccelerationMinimizerBlock(coreMatrix, splineDurations[splineId]);
    hessian_.block<SplineType::coefficientCount-2, SplineType::coefficientCount-2>(
        this->getSplineColumnIndex(splineId), this->getSplineColumnIndex(splineId)) =
            coreMatrix*weightMinAccel;

    // Regularization for position and velocity (cfMat must be positive definite).
    const unsigned int accelIdx = this->getCoeffIndex(splineId, SplineType::coefficientCount-2);
    hessian_.block<2,2>(accelIdx, accelIdx) += 1e-7*Eigen::Matrix2d::Identity();
  }

  return success;
}

template <int splineOrder_>
bool PolynomialSplineContainerOpt<splineOrder_>::setUpOptimizationMatrices(unsigned int num_coeffs) {
  // Inequality constraints.
  inequalityConstraintJacobian_.setZero(0, num_coeffs);
  inequalityConstraintMinValues_.setZero(0);

  // Add quadratic problem.
  costFunction_->setGlobalHessian(hessian_.sparseView());
  costFunction_->setLinearTerm(linearTerm_);
  functionConstraints_->setGlobalEqualityConstraintJacobian(this->equalityConstraintJacobian_.sparseView());
  functionConstraints_->setEqualityConstraintTargetValues(this->equalityConstraintTargetValues_);
  functionConstraints_->setGlobalInequalityConstraintJacobian(inequalityConstraintJacobian_.sparseView());
  functionConstraints_->setInequalityConstraintMinValues(inequalityConstraintMinValues_);
  functionConstraints_->setInequalityConstraintMaxValues(inequalityConstraintMinValues_);

  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainerOpt<splineOrder_>::solveQP(Eigen::VectorXd& coeffs) {
  bool success = true;

  double cost = 0.0;
  numopt_common::ParameterizationIdentity params(coeffs);
  success &= minimizer_->minimize(quadraticProblem_.get(), params, cost);

  if (!success) {
    std::cout << "[PolynomialSplineContainerOpt::solveQP] Failed to solve optimization!" << std::endl;
  } else {
    coeffs = params.getParams();
  }

  return success;
}


} /* namespace */
