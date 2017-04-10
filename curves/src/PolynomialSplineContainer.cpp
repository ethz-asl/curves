/*
 * PolynomialSplineContainer.cpp
 *
 *  Created on: Dec 8, 2014
 *      Author: C. Dario Bellicoso, Peter Fankhauser
 */

#include "curves/PolynomialSplineContainer.hpp"

// std
#include <iostream>

// boost
#include <boost/math/special_functions/pow.hpp>

namespace curves {

constexpr double PolynomialSplineContainer::undefinedValue;

// Import spline types from class.
using SplineType = PolynomialSplineContainer::SplineType;
using SplineList = PolynomialSplineContainer::SplineList;

PolynomialSplineContainer::PolynomialSplineContainer():
    timeOffset_(0.0),
    containerTime_(0.0),
    containerDuration_(0.0),
    activeSplineIdx_(0)
{
  // Make sure that the container is correctly emptied
  reset();
}


PolynomialSplineContainer::~PolynomialSplineContainer()
{
  // TODO Auto-generated destructor stub
}


bool PolynomialSplineContainer::advance(double dt)
{
  // Check container size
//  if (containerSize == 0 || containerSize == activeSplineIdx_) {
////    throw std::runtime_error("splinecontainer::advance");
//    return false;
//  }

  if (splines_.empty() || containerTime_ >= containerDuration_ || activeSplineIdx_ == splines_.size()) {
    return false;
  }

//  // Advance time from 0 to tf
//  if (time_ < splines_[activeSplineIdx_].getSplineDuration()) {
//    time_ += dt;
//  } else {
//    // Reset time
//    time_ = 0.0;
//    timeOffset_ += splines_[activeSplineIdx_].getSplineDuration();
//    activeSplineIdx_++;
//  }

  containerTime_ += dt;

  if ((containerTime_ - timeOffset_ >= splines_[activeSplineIdx_].getSplineDuration())) {
    if (activeSplineIdx_ < (splines_.size() - 1)) {
      timeOffset_ += splines_[activeSplineIdx_].getSplineDuration();
    }
    activeSplineIdx_++;
  }

  return true;
}

void PolynomialSplineContainer::setContainerTime(double t)
{
  containerTime_ = t;
  double timeOffset;
  activeSplineIdx_ = getActiveSplineIndexAtTime(t, timeOffset);
}

/*
 * aijh:
 *  i --> spline id (1,...,n)
 *  j --> spline coefficient aj (a5,...,a1,a0)
 *  h --> dimX, dimY
 *
 * Coefficient vector is:
 *    q = [a15x a14x ... a10x a15y ... a10y a25x ... a20y ... an5x ... an0y]
 */
int PolynomialSplineContainer::getCoeffIndex(int splineIdx, int aIdx) const {
  return splineIdx*6 + aIdx;
}

int PolynomialSplineContainer::getSplineColumnIndex(int splineIdx) const
{
  return getCoeffIndex(splineIdx, 0);
}

void PolynomialSplineContainer::setData(const std::vector<double>& knotPositions,
                                        const std::vector<double>& knotValues,
                                        double initialVelocity, double initialAcceleration,
                                        double finalVelocity, double finalAcceleration) {
  reset();

  const unsigned int num_splines = knotPositions.size()-1;
  constexpr auto num_coeffs_spline = SplineType::coefficientCount;
  const unsigned int num_coeffs = num_splines*num_coeffs_spline;
  const unsigned int num_knots = knotPositions.size();

  const unsigned int num_initial_constraints = 3;
  const unsigned int num_final_constraints = 3;

  const unsigned int num_constraints = (num_splines-1)*4 + num_initial_constraints + num_final_constraints;

  std::vector<double> tfs;// (num_splines);
  for (unsigned int i=0; i<num_splines; i++) {
    tfs.push_back(knotPositions[i+1]-knotPositions[i]);
  }

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_constraints, num_coeffs);
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(num_coeffs);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(num_constraints);

  // time containers
  SplineType::EigenTimeVectorType timeVec, dTimeVec, ddTimeVec;
  SplineType::EigenTimeVectorType timeVecTf, dTimeVecTf, ddTimeVecTf;

  SplineType::getTimeVector(timeVec, 0.0);
  SplineType::getdTimeVector(dTimeVec, 0.0);
  SplineType::getddTimeVector(ddTimeVec, 0.0);

  // define convenience quantities
//  int dimX = 0;
//  int a5 = 0;
//  int physPos = 0;
//  int physVel = 1;
//  int physAcc = 2;


  int constraintIdx = 0;
//  A.block(0, getCoeffIndex(1, a5),1,6) = timeVec;    // a_10_x
  A.block(constraintIdx, getSplineColumnIndex(0), 1, num_coeffs_spline) = timeVec;
  b(constraintIdx) = knotValues[0];
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(0), 1, num_coeffs_spline) = dTimeVec;   // a_11_x
  b(constraintIdx) = initialVelocity; // initial velocity
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(0), 1, num_coeffs_spline) = ddTimeVec;  // a_12_x
  b(constraintIdx) = initialAcceleration; // initial acceleration
  constraintIdx++;

  // Final conditions
//  int lastSplineId = num_splines; // doesnt have to be -1
  const double tf = tfs.back();
//  int rows = A.rows();

  SplineType::getTimeVector(timeVec, tf);
  SplineType::getdTimeVector(dTimeVec, tf);
  SplineType::getddTimeVector(ddTimeVec, tf);

  A.block(constraintIdx, getSplineColumnIndex(num_splines-1), 1, num_coeffs_spline) = timeVec;
  b(constraintIdx) = knotValues.back();
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(num_splines-1), 1, num_coeffs_spline) = dTimeVec;
  b(constraintIdx) = finalVelocity;
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(num_splines-1), 1, num_coeffs_spline) = ddTimeVec;
  b(constraintIdx) = finalAcceleration;
  constraintIdx++;
  /***************************/


  /**********************************
   * Set spline junction conditions *
   **********************************/
  for (size_t k=0; k<num_splines-1; k++) {

    const int prevSplineId = k;
    const int nextSplineId = k+1;

    const double tf = tfs[k];

    SplineType::getTimeVector(timeVec, 0.0);
    SplineType::getdTimeVector(dTimeVec, 0.0);
    SplineType::getddTimeVector(ddTimeVec, 0.0);

    SplineType::getTimeVector(timeVecTf, tf);
    SplineType::getdTimeVector(dTimeVecTf, tf);
    SplineType::getddTimeVector(ddTimeVecTf, tf);

    A.block(constraintIdx, getSplineColumnIndex(prevSplineId), 1, num_coeffs_spline) = timeVecTf;
    b(constraintIdx) = knotValues[k+1];
    constraintIdx++;

    A.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, num_coeffs_spline) = timeVec;
    b(constraintIdx) = knotValues[k+1];
    constraintIdx++;

    A.block(constraintIdx, getSplineColumnIndex(prevSplineId), 1, num_coeffs_spline) = dTimeVecTf;
    A.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, num_coeffs_spline) = -dTimeVec;
    b(constraintIdx) = 0.0;
    constraintIdx++;

    A.block(constraintIdx, getSplineColumnIndex(prevSplineId), 1, num_coeffs_spline) = ddTimeVecTf;
    A.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, num_coeffs_spline) = -ddTimeVec;
    b(constraintIdx) = 0.0;
    constraintIdx++;
  }
  /**********************************/

  coeffs = A.colPivHouseholderQr().solve(b);

  SplineType spline;
  SplineType::SplineCoefficients coefficients;

  for (unsigned int i = 0; i <num_splines; i++) {
    Eigen::Map<Eigen::VectorXd>(coefficients.data(), num_coeffs_spline, 1) = coeffs.segment<num_coeffs_spline>(getSplineColumnIndex(i));
    spline.setCoefficientsAndDuration(coefficients, tfs[i]);
    this->addSpline(spline);
  }

}

int PolynomialSplineContainer::getActiveSplineIndex() const
{
  return activeSplineIdx_;
}

bool PolynomialSplineContainer::addSpline(const SplineType& spline)
{
  splines_.push_back(spline);
  containerDuration_ += spline.getSplineDuration();
  return true;
}

bool PolynomialSplineContainer::reset()
{
  splines_.clear();
  activeSplineIdx_ = 0;
  containerDuration_ = 0.0;
  resetTime();
  return true;
}

bool PolynomialSplineContainer::resetTime()
{
  timeOffset_ = 0.0;
  containerTime_ = 0.0;
  activeSplineIdx_ = 0;
  return true;
}

double PolynomialSplineContainer::getContainerDuration() const
{
  return containerDuration_;
}

SplineType* PolynomialSplineContainer::getSpline(int splineIndex)
{
  return &splines_.at(splineIndex);
}

double PolynomialSplineContainer::getContainerTime() const
{
  return containerTime_;
}

bool PolynomialSplineContainer::isEmpty() const
{
  return splines_.empty();
}

double PolynomialSplineContainer::getPosition() const
{
  if (splines_.empty()) return 0.0;
  if (activeSplineIdx_ == splines_.size())
    return splines_.at(activeSplineIdx_ - 1).getPositionAtTime(containerTime_ - timeOffset_);
  return splines_.at(activeSplineIdx_).getPositionAtTime(containerTime_ - timeOffset_);
}


double PolynomialSplineContainer::getVelocity() const {
  if (splines_.empty()) return 0.0;
  if (activeSplineIdx_ == splines_.size())
    return splines_.at(activeSplineIdx_ - 1).getVelocityAtTime(containerTime_ - timeOffset_);
  return splines_.at(activeSplineIdx_).getVelocityAtTime(containerTime_ - timeOffset_);
}


double PolynomialSplineContainer::getAcceleration() const
{
  if (splines_.empty()) return 0.0;
  if (activeSplineIdx_ == splines_.size())
    return splines_.at(activeSplineIdx_ - 1).getAccelerationAtTime(containerTime_ - timeOffset_);
  return splines_.at(activeSplineIdx_).getAccelerationAtTime(containerTime_ - timeOffset_);
}

double PolynomialSplineContainer::getPositionAtTime(double t) const
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

int PolynomialSplineContainer::getActiveSplineIndexAtTime(double t, double& timeOffset) const {
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

double PolynomialSplineContainer::getVelocityAtTime(double t) const
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


double PolynomialSplineContainer::getAccelerationAtTime(double t) const
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

double PolynomialSplineContainer::getEndPosition() const
{
  double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
  return splines_.at(splines_.size() - 1).getPositionAtTime(lastSplineDuration);
}

double PolynomialSplineContainer::getEndVelocity() const
{
  double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
  return splines_.at(splines_.size() - 1).getVelocityAtTime(lastSplineDuration);
}

double PolynomialSplineContainer::getEndAcceleration() const
{
  double lastSplineDuration = splines_.at(splines_.size() - 1).getSplineDuration();
  return splines_.at(splines_.size() - 1).getAccelerationAtTime(lastSplineDuration);
}

const SplineList& PolynomialSplineContainer::getSplines() const {
  return splines_;
}

} /* namespace */
