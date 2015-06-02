/*
 * SlerpSE3Curve.cpp
 *
 *  Created on: Oct 10, 2014
 *      Author: Renaud Dube, Abel Gawel, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <curves/SlerpSE3Curve.hpp>
#include <iostream>

namespace curves {

SlerpSE3Curve::SlerpSE3Curve() : SE3Curve() {}

SlerpSE3Curve::~SlerpSE3Curve() {}

void SlerpSE3Curve::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=========== Slerp SE3 CURVE =============" << std::endl;
  std::cout << str << std::endl;
  std::cout << "num of coefficients: " << manager_.size() << std::endl;
  std::cout << "dimension: " << 6 << std::endl;
  std::stringstream ss;
  std::vector<Key> keys;
  std::vector<Time> times;
  manager_.getTimes(&times);
  manager_.getKeys(&keys);
  std::cout << "curve defined between times: " << manager_.getMinTime() <<
      " and " << manager_.getMaxTime() <<std::endl;
  std::cout <<"=========================================" <<std::endl;
  for (size_t i = 0; i < manager_.size(); i++) {
    ss << "coefficient " << keys[i] << ": ";
//    gtsam::traits<Coefficient>::Print(manager_.getCoefficientByKey(keys[i]),ss.str());
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
}

Time SlerpSE3Curve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time SlerpSE3Curve::getMinTime() const {
  return manager_.getMinTime();
}

bool SlerpSE3Curve::isEmpty() const {
  std::vector<Time> outTimes;
  manager_.getTimes(&outTimes);
  return outTimes.empty();
}

int SlerpSE3Curve::size() const {
  return manager_.size();
}

void SlerpSE3Curve::fitCurve(const std::vector<Time>& times,
                             const std::vector<ValueType>& values,
                             std::vector<Key>* outKeys) {
  assert(times.size() == values.size());
  if(times.size() > 0) {
    manager_.insertCoefficients(times, values, outKeys);
  }
}

void SlerpSE3Curve::extend(const std::vector<Time>& times,
                           const std::vector<ValueType>& values,
                           std::vector<Key>* outKeys) {

  if (times.size() != values.size())
    std::cerr << "number of times and number of coefficients don't match" << std::endl;

  manager_.insertCoefficients(times, values, outKeys);
}

typename SlerpSE3Curve::DerivativeType
SlerpSE3Curve::evaluateDerivative(
    Time time, unsigned derivativeOrder) const
{
  std::cerr << "Not implemented." << std::endl;

  // time is out of bound --> error
  if(time < this->getMinTime()) std::cerr << "Time out of bounds" << std::endl;
  if(time > this->getMaxTime()) std::cerr << "Time out of bounds" << std::endl;

  Eigen::VectorXd dCoeff;
  Time dt;
  CoefficientIter rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  if(!success) std::cerr << "Unable to get the coefficients at time " << time << std::endl;
  // first derivative
  // TODO!
//  if (derivativeOrder == 1) {
//    //todo Verify this
//    dCoeff = gtsam::traits<Coefficient>::Local(rval1->second.coefficient,rval0->second.coefficient);
//    dt = rval1->first - rval0->first;
//    return dCoeff/dt;
//    // order of derivative > 1 returns vector of zeros
//  } else {
//    const int dimension = gtsam::traits<Coefficient>::dimension;
//    return Eigen::VectorXd::Zero(dimension,1);
//  }
}

/// \brief \f[T^{\alpha}\f]
SE3 transformationPower(SE3 T, double alpha)
{
  SO3 R(T.getRotation());
  SE3::Position t(T.getPosition());

  AngleAxis angleAxis(R);
  angleAxis.setUnique();
  angleAxis.setAngle(angleAxis.angle() * alpha);
  angleAxis.setUnique();

  t *= alpha;
  R = angleAxis;
  return SE3(t, R);
}

/// \brief \f[A*B\f]
SE3 composeTransformations(SE3 A, SE3 B)
{
  SE3::Position position = A.transform(B.getPosition());
  SE3::Rotation rotation = A.getRotation() * B.getRotation();
  SE3 transformation(position, rotation);
  return transformation;
}

/// \brief \f[T^{-1}\f]
SE3 inverseTransformation(SE3 T)
{
  SE3::Position position = -T.getPosition();
  SE3::Rotation rotation = T.getRotation().inverted();
  SE3 inverse(position, rotation);
  return inverse;
}

SE3 invertAndComposeImplementation(SE3 A, SE3 B)
{
  SE3 result = composeTransformations(inverseTransformation(A), B);
  return result;
}

SE3 SlerpSE3Curve::evaluate(Time time) const
{
  std::cerr << "Not implemented." << std::endl;

  // Check if the curve is only defined at this one time
  if (manager_.getMaxTime() == time && manager_.getMinTime() == time) {
    return manager_.coefficientBegin()->second.coefficient;
  } else {
    CoefficientIter a, b;
    bool success = manager_.getCoefficientsAt(time, &a, &b);
    if(!success) std::cerr << "Unable to get the coefficients at time " << time << std::endl;
    SE3 T_W_A = a->second.coefficient;
    SE3 T_W_B = b->second.coefficient;
    double alpha = double(time - a->first)/double(b->first - a->first);

    // TODO.
//    //Implementation of T_W_I = T_W_A*exp(alpha*log(inv(T_W_A)*T_W_B))
//    using namespace kindr;
//    SE3 T_A_B = composeTransformations(inverseTransformation(T_W_A), T_W_B);
//    T_A_B.logarithmicMap();
//    gtsam::Vector6 log_T_A_B = transformationLogImplementation(T_A_B, boost::none);
//    gtsam::Vector6 log_T_A_I = vectorScalingImplementation<int(6)>(log_T_A_B, alpha, boost::none, boost::none);
//    SE3 T_A_I = transformationExpImplementation(log_T_A_I, boost::none);
//    return composeImplementation(T_W_A, T_A_I, boost::none, boost::none);
  }
}

void SlerpSE3Curve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  std::cerr << "Not implemented." << std::endl;
}

/// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityA(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityB(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityA(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityB(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistA(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistB(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeA(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeB(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}

} // namespace
