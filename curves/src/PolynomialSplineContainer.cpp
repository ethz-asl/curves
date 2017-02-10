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
 * Time vector tau(tk) is defined as:
 *  tau(tk) = [ tk^5  tk^4  tk^3  tk^2  tk  1].'
 */
void getTimeVector(Eigen::Matrix<double, 1, 6>& timeVec, double t_k)
{
  timeVec(5) = 1.0;
  timeVec(4) = t_k;
  timeVec(3) = t_k * timeVec(4);
  timeVec(2) = t_k * timeVec(3);
  timeVec(1) = t_k * timeVec(2);
  timeVec(0) = t_k * timeVec(1);
}

/*
 * Time vector dtau(tk) is defined as:
 *  dtau(tk) = [ 5tk^4  4tk^3  3tk^2  2tk  1  0].'
 */
void getdTimeVector(Eigen::Matrix<double, 1, 6>& timeVec, double t_k)
{
//  timeVec(0) = 5.0*boost::math::pow<4>(t_k);
//  timeVec(1) = 4.0*boost::math::pow<3>(t_k);
//  timeVec(2) = 3.0*boost::math::pow<2>(t_k);
  timeVec(0) = 5.0 * t_k * t_k * t_k * t_k;
  timeVec(1) = 4.0 * t_k * t_k * t_k;
  timeVec(2) = 3.0 * t_k * t_k;
  timeVec(3) = 2.0 * t_k;
  timeVec(4) = 1.0;
  timeVec(5) = 0.0;
}

/*
 * Time vector ddtau(tk) is defined as:
 *  ddtau(tk) = [ 20tk^3  12tk^2  6tk  2  0  0].'
 */
void getddTimeVector(Eigen::Matrix<double, 1, 6>& timeVec, double t_k)
{
//  timeVec(0) = 20.0*boost::math::pow<3>(t_k);
//  timeVec(1) = 12.0*boost::math::pow<2>(t_k);
  timeVec(0) = 20.0 * t_k * t_k * t_k;
  timeVec(1) = 12.0 * t_k * t_k;
  timeVec(2) = 6.0 * t_k;
  timeVec(3) = 2.0;
  timeVec(4) = 0.0;
  timeVec(5) = 0.0;
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
int getCoeffIndex(int splineIdx, int aIdx)
{
  int splineOffset = (splineIdx - 1) * 6;
  int idx = splineOffset + aIdx;

  return idx;
}

inline int getSplineColumnIndex(int splineIdx)
{
  return getCoeffIndex(splineIdx, 0);
}

void PolynomialSplineContainer::setData(const std::vector<double>& knotPositions,
                                        const std::vector<double>& knotValues,
                                        double initialVelocity, double initialAcceleration,
                                        double finalVelocity, double finalAcceleration)
{
//  for (int k =0; k< knotPositions.size(); k++) {
//    std::cout << "pos: " << knotPositions[k] << " val: " << knotValues[k] << std::endl;
//  }
  reset();

  unsigned int num_splines = knotPositions.size()-1;
  unsigned int num_coeffs_spline = 6;
  unsigned int num_coeffs = num_splines*num_coeffs_spline;
  unsigned int num_knots = knotPositions.size();
//  int num_constraints = knotPositions.size()*3+(knotPositions.size()-1);

  unsigned int num_initial_constraints = 3;
  unsigned int num_final_constraints = 3;

  unsigned int num_constraints = (num_knots-2)*4 + num_initial_constraints + num_final_constraints;

  std::vector<double> tfs;// (num_splines);
  for (unsigned int i=0; i<num_splines; i++) {
    tfs.push_back(knotPositions[i+1]-knotPositions[i]);
//    tfs[k] = knotPositions[k+1]-knotPositions[k];
//    std::cout << "duration: " << tfs[i] << " knot: " << i << std::endl;
  }

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_constraints, num_coeffs);
  Eigen::VectorXd coeffs = Eigen::VectorXd::Zero(num_coeffs);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(num_constraints);

  // time containers
  Eigen::Matrix<double,1,6> timeVec, dTimeVec, ddTimeVec;
  Eigen::Matrix<double,1,6> timeVecTf, dTimeVecTf, ddTimeVecTf;

  getTimeVector(timeVec, 0.0);
  getdTimeVector(dTimeVec, 0.0);
  getddTimeVector(ddTimeVec, 0.0);

  // define convenience quantities
//  int dimX = 0;
//  int a5 = 0;
//  int physPos = 0;
//  int physVel = 1;
//  int physAcc = 2;


  int constraintIdx = 0;
//  A.block(0, getCoeffIndex(1, a5),1,6) = timeVec;    // a_10_x
  A.block(constraintIdx, getSplineColumnIndex(1), 1, num_coeffs_spline) = timeVec;
  b(constraintIdx) = knotValues[0];
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(1), 1, num_coeffs_spline) = dTimeVec;   // a_11_x
  b(constraintIdx) = initialVelocity; // initial velocity
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(1), 1, num_coeffs_spline) = ddTimeVec;  // a_12_x
  b(constraintIdx) = initialAcceleration; // initial acceleration
  constraintIdx++;

  // Final conditions
//  int lastSplineId = num_splines; // doesnt have to be -1
  double tf = tfs.back();
//  int rows = A.rows();

  getTimeVector(timeVec, tf);
  getdTimeVector(dTimeVec, tf);
  getddTimeVector(ddTimeVec, tf);

  A.block(constraintIdx, getSplineColumnIndex(num_splines), 1, num_coeffs_spline) = timeVec;
  b(constraintIdx) = knotValues.back();
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(num_splines), 1, num_coeffs_spline) = dTimeVec;
  b(constraintIdx) = finalVelocity;
  constraintIdx++;

  A.block(constraintIdx, getSplineColumnIndex(num_splines), 1, num_coeffs_spline) = ddTimeVec;
  b(constraintIdx) = finalAcceleration;
  constraintIdx++;
  /***************************/


  /**********************************
   * Set spline junction conditions *
   **********************************/
  for (size_t k=1; k<=num_knots-2; k++) {

    int prevSplineId = k;
    int nextSplineId = k+1;

    double tf = tfs[k-1];

    getTimeVector(timeVec, 0.0);
    getdTimeVector(dTimeVec, 0.0);
    getddTimeVector(ddTimeVec, 0.0);

    getTimeVector(timeVecTf, tf);
    getdTimeVector(dTimeVecTf, tf);
    getddTimeVector(ddTimeVecTf, tf);

    A.block(constraintIdx, getSplineColumnIndex(prevSplineId), 1, num_coeffs_spline) = timeVecTf;
    b(constraintIdx) = knotValues[k];
//    std::cout << "knot val: " << knotValues[k] << std::endl;
    constraintIdx++;

    A.block(constraintIdx, getSplineColumnIndex(nextSplineId), 1, num_coeffs_spline) = timeVec;
    b(constraintIdx) = knotValues[k];
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


  PolynomialSplineQuintic spline;
  std::vector<double> coefficients;

//  std::cout << "number of splines: " << num_splines << std::endl;
  for (unsigned int i = 0; i <num_splines; i++) {
    coefficients.clear();
    for (int k = num_coeffs_spline-1; k >= 0; k--) {
      coefficients.push_back( static_cast<double>(coeffs( getSplineColumnIndex(i+1)+k ) ));
    }
    spline.setCoeffsAndDuration(coefficients, tfs[i]);
    this->addSpline(spline);
  }

//  Eigen::IOFormat CleanFmt(2, 0, ",","\n", "[", "]");
//  std::cout  << "A:\n"  << A.format(CleanFmt) << std::endl;
//  std::cout  << "b:\n"  << b.format(CleanFmt) << std::endl;
//  std::cout << "coeffs: " << coeffs << std::endl;

}

int PolynomialSplineContainer::getActiveSplineIndex() const
{
  return activeSplineIdx_;
}

bool PolynomialSplineContainer::addSpline(const PolynomialSplineQuintic& spline)
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

PolynomialSplineBase* PolynomialSplineContainer::getSpline(int splineIndex)
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
//  std::cout << "splineIdx: " << activeSplineIdx_ << std::endl;
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


//double PolynomialSplineContainer::getVelocityAtTime(double t) {
//  throw std::runtime_error("not yet implemented");
//}

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

//double PolynomialSplineContainer::getAccelerationAtTime(double t) {
//  throw std::runtime_error("not yet implemented");
//}

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

const std::vector<PolynomialSplineQuintic>& PolynomialSplineContainer::getSplines() const {
  return splines_;
}

} /* namespace */
