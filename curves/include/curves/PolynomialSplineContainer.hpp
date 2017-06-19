/*
 * PolynomialSplineContainer.hpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#pragma once

#include "curves/polynomial_splines.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits>

namespace curves {

template <int splineOrder_>
class PolynomialSplineContainer {
 public:

  using SplineType = PolynomialSpline<splineOrder_>;
  using SplineList = std::vector<SplineType>;

  PolynomialSplineContainer():
      timeOffset_(0.0),
      containerTime_(0.0),
      containerDuration_(0.0),
      activeSplineIdx_(0),
      hessian_(),
      linearTerm_()
  {
    // Make sure that the container is correctly emptied
    reset();
  }


  ~PolynomialSplineContainer()
  {
    // TODO Auto-generated destructor stub
  }



  bool advance(double dt)
  {
    if (splines_.empty() || containerTime_ >= containerDuration_ || activeSplineIdx_ == splines_.size()) {
      return false;
    }

    containerTime_ += dt;

    // check if spline index needs to be increased
    if ((containerTime_ - timeOffset_ >= splines_[activeSplineIdx_].getSplineDuration())) {
      if (activeSplineIdx_ < (splines_.size() - 1)) {
        timeOffset_ += splines_[activeSplineIdx_].getSplineDuration();
      }
      activeSplineIdx_++;
    }

    return true;
  }


  void setContainerTime(double t)
  {
    containerTime_ = t;
    double timeOffset;
    activeSplineIdx_ = getActiveSplineIndexAtTime(t, timeOffset);
  }


  int getActiveSplineIndex() const
  {
    return activeSplineIdx_;
  }


  bool addSpline(const SplineType& spline)
  {
    splines_.push_back(spline);
    containerDuration_ += spline.getSplineDuration();
    return true;
  }


  bool addSpline(SplineType&& spline) {
    containerDuration_ += spline.getSplineDuration();
    splines_.emplace_back(spline);
    return true;
  }


  bool reset()
  {
    splines_.clear();
    activeSplineIdx_ = 0;
    containerDuration_ = 0.0;
    resetTime();
    return true;
  }


  bool resetTime()
  {
    timeOffset_ = 0.0;
    containerTime_ = 0.0;
    activeSplineIdx_ = 0;
    return true;
  }


  double getContainerDuration() const
  {
    return containerDuration_;
  }


  SplineType* getSpline(int splineIndex)
  {
    return &splines_.at(splineIndex);
  }


  double getContainerTime() const
  {
    return containerTime_;
  }


  bool isEmpty() const
  {
    return splines_.empty();
  }


  double getPosition() const
  {
    if (splines_.empty()) return 0.0;
    if (activeSplineIdx_ == splines_.size())
      return splines_.at(activeSplineIdx_ - 1).getPositionAtTime(containerTime_ - timeOffset_);
    return splines_.at(activeSplineIdx_).getPositionAtTime(containerTime_ - timeOffset_);
  }



  double getVelocity() const {
    if (splines_.empty()) return 0.0;
    if (activeSplineIdx_ == splines_.size())
      return splines_.at(activeSplineIdx_ - 1).getVelocityAtTime(containerTime_ - timeOffset_);
    return splines_.at(activeSplineIdx_).getVelocityAtTime(containerTime_ - timeOffset_);
  }



  double getAcceleration() const
  {
    if (splines_.empty()) return 0.0;
    if (activeSplineIdx_ == splines_.size())
      return splines_.at(activeSplineIdx_ - 1).getAccelerationAtTime(containerTime_ - timeOffset_);
    return splines_.at(activeSplineIdx_).getAccelerationAtTime(containerTime_ - timeOffset_);
  }


  double getPositionAtTime(double t) const
  {
    double timeOffset = 0.0;
    int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);

    if (activeSplineIdx < 0) {
      return splines_.at(0).getPositionAtTime(0.0);
    }
    if (activeSplineIdx == splines_.size()) {
      return splines_.at(activeSplineIdx - 1).getPositionAtTime(
          splines_.at(activeSplineIdx - 1).getSplineDuration());
    }
    return splines_.at(activeSplineIdx).getPositionAtTime(t - timeOffset);
  }


  int getActiveSplineIndexAtTime(double t, double& timeOffset) const {
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


  double getVelocityAtTime(double t) const
  {
    double timeOffset = 0.0;
    int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);
    if (activeSplineIdx < 0) {
      return splines_.at(0).getVelocityAtTime(0.0);
    }
    if (activeSplineIdx == splines_.size()) {
      return splines_.at(activeSplineIdx - 1).getVelocityAtTime(
          splines_.at(activeSplineIdx - 1).getSplineDuration());
    }
    return splines_.at(activeSplineIdx).getVelocityAtTime(t - timeOffset);
  }



  double getAccelerationAtTime(double t) const
  {
    double timeOffset = 0.0;
    int activeSplineIdx = getActiveSplineIndexAtTime(t, timeOffset);
    if (activeSplineIdx < 0) {
      return splines_.at(0).getAccelerationAtTime(0.0);
    }
    if (activeSplineIdx == splines_.size()) {
      return splines_.at(activeSplineIdx - 1).getAccelerationAtTime(
          splines_.at(activeSplineIdx - 1).getSplineDuration());
    }
    return splines_.at(activeSplineIdx).getAccelerationAtTime(t - timeOffset);
  }


  double getEndPosition() const
  {
    double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
    return splines_.at(splines_.size() - 1).getPositionAtTime(lastSplineDuration);
  }


  double getEndVelocity() const
  {
    double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
    return splines_.at(splines_.size() - 1).getVelocityAtTime(lastSplineDuration);
  }


  double getEndAcceleration() const
  {
    double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
    return splines_.at(splines_.size() - 1).getAccelerationAtTime(lastSplineDuration);
  }


  const SplineList& getSplines() const {
    return splines_;
  }




  /*
   * !Minimize spline coefficients s.t. position, velocity and acceleration constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the second derivative).
   * If the spline order is smaller than 5, some boundary constraints are dropped.
   * Notice that the number of junction constraints (and thereby the degree of smoothness at the transition
   * between two adjacent splines) decreases with decreasing spline order.
   */
  bool setData(const std::vector<double>& knotDurations,
               const std::vector<double>& knotPositions,
               double initialVelocity,
               double initialAcceleration,
               double finalVelocity,
               double finalAcceleration) {

    bool success = true;

    success &= reset();

    // Set up optimization parameters
    const unsigned int num_splines = knotDurations.size()-1;
    constexpr auto num_coeffs_spline = SplineType::coefficientCount;
    const unsigned int num_coeffs = num_splines*num_coeffs_spline;
    const unsigned int num_junctions = num_splines-1;

    if (num_splines<1) {
      return false;
    }

    // total number of constraints
    constexpr unsigned int num_initial_constraints = 3;   // pos, vel, accel
    constexpr unsigned int num_final_constraints = 3;     // pos, vel, accel
    constexpr unsigned int num_constraint_junction = 4;   // pos (2x), vel, accel
    const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
    const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

    // dop constraints if necessary
    if (num_constraints>num_coeffs) {
      std::cout << "[PolynomialSplineContainer::setData] Number of equality constraints is larger than number of coefficients. Drop acceleration constraints!\n";
      return setData(knotDurations, knotPositions, initialVelocity, finalVelocity);
    }

    // Vector containing durations of splines
    std::vector<double> splineDurations(num_splines);
    for (unsigned int splineId=0; splineId<num_splines; splineId++) {
      splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

      if (splineDurations[splineId]<=0.0) {
        std::cout << "[PolynomialSplineContainer::setData] Invalid spline duration!\n";
        return false;
      }
    }

    // Initialize optimization matrices
    hessian_.setZero(num_constraints, num_coeffs);
    linearTerm_.setZero(num_constraints);
    Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(num_coeffs);
    unsigned int constraintIdx = 0;

    // Initial conditions
    Eigen::VectorXd initialConditions(num_initial_constraints);
    initialConditions(0) = knotPositions.front();
    initialConditions(1) = initialVelocity;
    initialConditions(2) = initialAcceleration;
    addInitialConditions(initialConditions, constraintIdx);

    // Final conditions
    Eigen::VectorXd finalConditions(num_final_constraints);
    finalConditions(0) = knotPositions.back();
    finalConditions(1) = finalVelocity;
    finalConditions(2) = finalAcceleration;
    addFinalConditions(finalConditions, constraintIdx,  splineDurations.back(), num_junctions);

    // Junction conditions
    addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

    if (num_constraints!=constraintIdx) {
      std::cout << "[PolynomialSplineContainer::setData] Wrong number of equality constraints. Num of equality constraints =  "
          << num_constraints << ", counter = " << constraintIdx << std::endl;
      return false;
    }


    // Minimize spline coefficients
    coeffs = hessian_.colPivHouseholderQr().solve(linearTerm_);

    //std::cout << "wrong hessian_ = \n" << hessian_ << std::endl;
    std::cout << "wrong coeffs = \n" << coeffs << std::endl;

    // Extract spline coefficients and add splines
    SplineType spline;
    typename SplineType::SplineCoefficients coefficients;

    for (unsigned int i = 0; i <num_splines; i++) {
      Eigen::Map<Eigen::VectorXd>(coefficients.data(), num_coeffs_spline, 1) = coeffs.segment<num_coeffs_spline>(getSplineColumnIndex(i));
      spline.setCoefficientsAndDuration(coefficients, splineDurations[i]);
      this->addSpline(spline);
    }

    return success;

  }

  /*
   * !Minimize spline coefficients s.t. position and velocity constraints are satisfied
   * (i.e., s.t. the spline conjunction is smooth up the first derivative).
   * If the spline order is smaller than 3, some boundary constraints are dropped.
   * Notice that the number of junction constraints (and thereby the degree of smoothness at the transition
   * between two adjacent splines) decreases with decreasing spline order.
   */
  bool setData(const std::vector<double>& knotDurations,
               const std::vector<double>& knotPositions,
               double initialVelocity,
               double finalVelocity) {

    bool success = true;

    success &= reset();

    // Set up optimization parameters
    const unsigned int num_splines = knotDurations.size()-1;
    constexpr auto num_coeffs_spline = SplineType::coefficientCount;
    const unsigned int num_coeffs = num_splines*num_coeffs_spline;
    const unsigned int num_junctions = num_splines-1;

    if (num_splines<1) {
      return false;
    }

    // total number of constraints
    constexpr unsigned int num_initial_constraints = 2;   // pos, vel
    constexpr unsigned int num_final_constraints = 2;     // pos, vel
    constexpr unsigned int num_constraint_junction = 3;   // pos (2x), vel
    const unsigned int num_junction_constraints = num_junctions*num_constraint_junction;
    const unsigned int num_constraints = num_junction_constraints + num_initial_constraints + num_final_constraints;

    // dop constraints if necessary
    if (num_constraints>num_coeffs) {
      std::cout << "[PolynomialSplineContainer::setData] Number of equality constraints is larger than number of coefficients. Drop velocity constraints!\n";
      return setData(knotDurations, knotPositions);
    }

    // Vector containing durations of splines
    std::vector<double> splineDurations(num_splines);
    for (unsigned int splineId=0; splineId<num_splines; splineId++) {
      splineDurations[splineId] = knotDurations[splineId+1]-knotDurations[splineId];

      if (splineDurations[splineId]<=0.0) {
        std::cout << "[PolynomialSplineContainer::setData] Invalid spline duration!\n";
        return false;
      }
    }

    // Initialize optimization matrices
    hessian_.setZero(num_constraints, num_coeffs);
    linearTerm_.setZero(num_constraints);
    Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(num_coeffs);
    unsigned int constraintIdx = 0;

    // Initial conditions
    Eigen::VectorXd initialConditions(num_initial_constraints);
    initialConditions(0) = knotPositions.front();
    initialConditions(1) = initialVelocity;
    addInitialConditions(initialConditions, constraintIdx);

    // Final conditions
    Eigen::VectorXd finalConditions(num_final_constraints);
    finalConditions(0) = knotPositions.back();
    finalConditions(1) = finalVelocity;
    addFinalConditions(finalConditions, constraintIdx,  splineDurations.back(), num_junctions);

    // Junction conditions
    addJunctionsConditions(splineDurations, knotPositions, constraintIdx, num_junctions);

    if (num_constraints!=constraintIdx) {
      std::cout << "[PolynomialSplineContainer::setData] Wrong number of equality constraints. Num of equality constraints =  "
          << num_constraints << ", counter = " << constraintIdx << std::endl;
      return false;
    }

    // Minimize spline coefficients
    coeffs = hessian_.colPivHouseholderQr().solve(linearTerm_);

    // Extract spline coefficients and add splines
    SplineType spline;
    typename SplineType::SplineCoefficients coefficients;

    for (unsigned int splineId = 0; splineId <num_splines; splineId++) {
      Eigen::Map<Eigen::VectorXd>(coefficients.data(), num_coeffs_spline, 1) = coeffs.segment<num_coeffs_spline>(getSplineColumnIndex(splineId));
      spline.setCoefficientsAndDuration(coefficients, knotDurations[splineId]);
      success &= this->addSpline(spline);
    }

    return success;

  }


  /*
   * ! Find linear part of the spline coefficients (a0, a1) s.t. position constraints are satisfied.
   * If the spline order is larger than 1, the remaining spline coefficients are set to zero.
   */
  bool setData(const std::vector<double>& knotDurations,
               const std::vector<double>& knotPositions) {

    if (splineOrder_==0) {
      return false;
    }

    bool success = true;
    const int num_splines = knotDurations.size()-1;
    constexpr auto num_coeffs_spline = SplineType::coefficientCount;
    typename SplineType::SplineCoefficients coefficients;

    if (num_splines<1) {
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





  static constexpr double undefinedValue = std::numeric_limits<double>::quiet_NaN();

 protected:


  /*
   * aijh:
   *  i --> spline id (1,...,n)
   *  j --> spline coefficient aj (a5,...,a1,a0)
   *  h --> dimX, dimY
   *
   * Coefficient vector is:
   *    q = [a15x a14x ... a10x a15y ... a10y a25x ... a20y ... an5x ... an0y]
   */

  int getCoeffIndex(int splineIdx, int aIdx) const {
    return splineIdx*(splineOrder_+1) + aIdx;
  }


  int getSplineColumnIndex(int splineIdx) const
  {
    return getCoeffIndex(splineIdx, 0);
  }


  void addInitialConditions(const Eigen::VectorXd& initialConditions,
                            unsigned int& constraintIdx) {

    // time container
    typename SplineType::EigenTimeVectorType timeVec;

    // initial position
    if (initialConditions.size()>0) {
      SplineType::getTimeVector(timeVec, 0.0);
      hessian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
      linearTerm_(constraintIdx) = initialConditions(0);
      constraintIdx++;
    }

    // initial velocity
    if (initialConditions.size()>1) {
      SplineType::getdTimeVector(timeVec, 0.0);
      hessian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
      linearTerm_(constraintIdx) = initialConditions(1);
      constraintIdx++;
    }

    // initial acceleration
    if (initialConditions.size()>2) {
      SplineType::getddTimeVector(timeVec, 0.0);
      hessian_.block(constraintIdx, getSplineColumnIndex(0), 1, SplineType::coefficientCount) = timeVec;
      linearTerm_(constraintIdx) = initialConditions(2);
      constraintIdx++;
    }
  }


  void addFinalConditions(const Eigen::VectorXd& finalConditions,
                          unsigned int& constraintIdx,
                          double lastSplineDuration,
                          unsigned int lastSplineId) {

    // time container
    typename SplineType::EigenTimeVectorType timeVec0;


    // initial position
    if (finalConditions.size()>0) {
      SplineType::getTimeVector(timeVec0, lastSplineDuration);
      hessian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec0;
      linearTerm_(constraintIdx) = finalConditions(0);
      constraintIdx++;
    }

    // initial velocity
    if (finalConditions.size()>1) {
      SplineType::getdTimeVector(timeVec0, lastSplineDuration);
      hessian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec0;
      linearTerm_(constraintIdx) = finalConditions(1);
      constraintIdx++;
    }

    // initial acceleration
    if (finalConditions.size()>2) {
      SplineType::getddTimeVector(timeVec0, lastSplineDuration);
      hessian_.block(constraintIdx, getSplineColumnIndex(lastSplineId), 1, SplineType::coefficientCount) = timeVec0;
      linearTerm_(constraintIdx) = finalConditions(2);
      constraintIdx++;
    }
  }


  void addJunctionsConditions(const std::vector<double>& splineDurations,
                              const std::vector<double>& knotPositions,
                              unsigned int& constraintIdx,
                              unsigned int num_junctions) {

    // time containers
    typename SplineType::EigenTimeVectorType timeVec0, dTimeVec0, ddTimeVec0;
    typename SplineType::EigenTimeVectorType timeVecTf, dTimeVecTf, ddTimeVecTf;

    // get time container at zero time
    SplineType::getTimeVector(timeVec0, 0.0);
    SplineType::getdTimeVector(dTimeVec0, 0.0);
    SplineType::getddTimeVector(ddTimeVec0, 0.0);


    for (unsigned int splineId=0; splineId<num_junctions; splineId++) {
      const unsigned int nextSplineId = splineId+1;

      // get time container at spline duration
      SplineType::getTimeVector(timeVecTf, splineDurations[splineId]);
      SplineType::getdTimeVector(dTimeVecTf, splineDurations[splineId]);
      SplineType::getddTimeVector(ddTimeVecTf, splineDurations[splineId]);

      // smooth position transition with fixed positions
      hessian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  timeVecTf;
      linearTerm_(constraintIdx) = knotPositions[nextSplineId];
      constraintIdx++;

      hessian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) =  timeVec0;
      linearTerm_(constraintIdx) = knotPositions[nextSplineId];
      constraintIdx++;

      // smooth velocity transition
      hessian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  dTimeVecTf;
      hessian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) = -dTimeVec0;
      linearTerm_(constraintIdx) = 0.0;
      constraintIdx++;

      // smooth acceleration transition
      hessian_.block(constraintIdx, getSplineColumnIndex(splineId),     1, SplineType::coefficientCount) =  ddTimeVecTf;
      hessian_.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, SplineType::coefficientCount) = -ddTimeVec0;
      linearTerm_(constraintIdx) = 0.0;
      constraintIdx++;
    }
  }

  //! Conjunction of smoothly interconnected splines.
  SplineList splines_;

  //! Helper variable.
  double timeOffset_;

  //! Current time of spline conjunction.
  double containerTime_;

  //! Total duration of spline conjunction.
  double containerDuration_;

  //! Spline index currently active.
  int activeSplineIdx_;

  // Hessian of quadratic program.
  Eigen::MatrixXd hessian_;

  // Linear term of quadratic program.
  Eigen::VectorXd linearTerm_;
};

} /* namespace */
