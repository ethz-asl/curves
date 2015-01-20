/*
 * @file SlerpSE3Curve.cpp
 * @date Oct 10, 2014
 * @author Renaud Dube, Abel Gawel
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
    gtsam::traits<Coefficient>::Print(manager_.getCoefficientByKey(keys[i]),ss.str());
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
  CHECK_EQ(times.size(), values.size());
  if(times.size() > 0) {
    manager_.insertCoefficients(times,values, outKeys);
  }
}

void SlerpSE3Curve::extend(const std::vector<Time>& times,
                           const std::vector<ValueType>& values,
                           std::vector<Key>* outKeys) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";

  manager_.insertCoefficients(times, values, outKeys);
}

typename SlerpSE3Curve::DerivativeType
SlerpSE3Curve::evaluateDerivative(Time time,
                                  unsigned derivativeOrder) const {

  // time is out of bound --> error
  CHECK_GE(time, this->getMinTime()) << "Time out of bounds";
  CHECK_LE(time, this->getMaxTime()) << "Time out of bounds";

  Eigen::VectorXd dCoeff;
  Time dt;
  CoefficientIter rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  // first derivative
  if (derivativeOrder == 1) {
    //todo Verify this
    dCoeff = gtsam::traits<Coefficient>::Local(rval1->second.coefficient,rval0->second.coefficient);
    dt = rval1->first - rval0->first;
    return dCoeff/dt;
  } else { // order of derivative > 1 returns vector of zeros
    const int dimension = gtsam::traits<Coefficient>::dimension;
    return Eigen::VectorXd::Zero(dimension,1);
  }
}

/// \brief \f[T^{\alpha}\f]
SE3 transformationPower(SE3  T, double alpha) {
  SO3 R(T.getRotation());
  SE3::Position t(T.getPosition());

  AngleAxis angleAxis(R);
  angleAxis.setUnique();
  angleAxis.setAngle( angleAxis.angle()*alpha);
  angleAxis.setUnique();

  return SE3(SO3(angleAxis),(t*alpha).eval());
}

/// \brief \f[A*B\f]
SE3 composeTransformations(SE3 A, SE3 B) {
  return A*B;
}

/// \brief \f[T^{-1}\f]
SE3 inverseTransformation(SE3 T) {
  return T.inverted();
}

/// \brief forms slerp interpolation into a binary expression with 2 leafs and binds alpha into it,
///        uses break down of expression into its operations
///        \f[ T = A(A^{-1}B)^{\alpha} \f]
gtsam::Expression<typename SlerpSE3Curve::ValueType>
SlerpSE3Curve::getValueExpression(const Time& time) const {
  typedef typename SlerpSE3Curve::ValueType ValueType;
  using namespace gtsam;
  CoefficientIter rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  Expression<ValueType> leaf1(rval0->second.key);
  Expression<ValueType> leaf2(rval1->second.key);
  double alpha = double(time - rval0->first)/double(rval1->first - rval0->first);

  if (alpha == 0) {
    return leaf1;
  } else if (alpha == 1) {
    return leaf2;
  } else {
    return slerp(leaf1, leaf2, alpha);
  }
}

gtsam::Expression<typename SlerpSE3Curve::DerivativeType>
SlerpSE3Curve::getDerivativeExpression(const Time& time, unsigned derivativeOrder) const {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

SE3 SlerpSE3Curve::evaluate(Time time) const {
  CoefficientIter a, b;
  bool success = manager_.getCoefficientsAt(time, &a, &b);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  double alpha = double(time - a->first)/double(b->first - a->first);

  SE3 delta = composeTransformations(inverseTransformation(a->second.coefficient),b->second.coefficient);

  return composeTransformations(a->second.coefficient, transformationPower(delta,alpha));
}


void SlerpSE3Curve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

/// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}

void SlerpSE3Curve::initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) const {
  manager_.initializeGTSAMValues(keys, values);
}

void SlerpSE3Curve::initializeGTSAMValues(gtsam::Values* values) const {
  manager_.initializeGTSAMValues(values);
}

void SlerpSE3Curve::updateFromGTSAMValues(const gtsam::Values& values) {
  manager_.updateFromGTSAMValues(values);
}

} // namespace curves
