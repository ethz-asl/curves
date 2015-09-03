/*
 * @file SlerpSE3Curve.cpp
 * @date Oct 10, 2014
 * @author Renaud Dube, Abel Gawel
 */

#include <curves/SlerpSE3Curve.hpp>
#include <iostream>

#include "gtsam/nonlinear/ExpressionFactor.h"

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
  double sum_dp = 0;
  Eigen::Vector3d p1, p2;
  for(size_t i = 0; i < times.size()-1; ++i) {
    p1 = evaluate(times[i]).getPosition();
    p2 = evaluate(times[i+1]).getPosition();
    sum_dp += (p1-p2).norm();
  }
  std::cout << "average dt between coefficients: " << (manager_.getMaxTime() -manager_.getMinTime())  / (times.size()-1) << " ns." << std::endl;
  std::cout << "average distance between coefficients: " << sum_dp / double((times.size()-1))<< " m." << std::endl;
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
  return manager_.empty();
}

int SlerpSE3Curve::size() const {
  return manager_.size();
}

void SlerpSE3Curve::fitCurve(const std::vector<Time>& times,
                             const std::vector<ValueType>& values,
                             std::vector<Key>* outKeys) {
  CHECK_EQ(times.size(), values.size());
  if(times.size() > 0) {
    clear();
    manager_.insertCoefficients(times,values, outKeys);
  }
}

void SlerpSE3Curve::setCurve(const std::vector<Time>& times,
              const std::vector<ValueType>& values) {
  CHECK_EQ(times.size(), values.size());
  if(times.size() > 0) {
    manager_.insertCoefficients(times,values);
  }
}


void SlerpSE3Curve::extend(const std::vector<Time>& times,
                           const std::vector<ValueType>& values,
                           std::vector<Key>* outKeys) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  slerpPolicy_.extend<SlerpSE3Curve, ValueType>(times, values, this, outKeys);
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
    // order of derivative > 1 returns vector of zeros
  } else {
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
  // Check if the curve is only defined at this one time
  if (manager_.getMaxTime() == time && manager_.getMinTime() == time) {
    return manager_.coefficientBegin()->second.coefficient;
  } else {
    if (time == manager_.getMaxTime()) {
      // Efficient evaluation of a curve end
      return (--manager_.coefficientEnd())->second.coefficient;
    } else {
      CoefficientIter a, b;
      bool success = manager_.getCoefficientsAt(time, &a, &b);
      CHECK(success) << "Unable to get the coefficients at time " << time;
      SE3 T_W_A = a->second.coefficient;
      SE3 T_W_B = b->second.coefficient;
      double alpha = double(time - a->first)/double(b->first - a->first);

      //Implementation of T_W_I = T_W_A*exp(alpha*log(inv(T_W_A)*T_W_B))
      using namespace kindr::minimal;
      SE3 T_A_B = invertAndComposeImplementation(T_W_A, T_W_B, boost::none, boost::none);
      gtsam::Vector6 log_T_A_B = transformationLogImplementation(T_A_B, boost::none);
      gtsam::Vector6 log_T_A_I = vectorScalingImplementation<int(6)>(log_T_A_B, alpha, boost::none, boost::none);
      SE3 T_A_I = transformationExpImplementation(log_T_A_I, boost::none);
      return composeImplementation(T_W_A, T_A_I, boost::none, boost::none);
    }
  }
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

void SlerpSE3Curve::initializeGTSAMValues(gtsam::KeySet keys, gtsam::Values* values) const {
  manager_.initializeGTSAMValues(keys, values);
}

void SlerpSE3Curve::initializeGTSAMValues(gtsam::Values* values) const {
  manager_.initializeGTSAMValues(values);
}

void SlerpSE3Curve::updateFromGTSAMValues(const gtsam::Values& values) {
  manager_.updateFromGTSAMValues(values);
}

void SlerpSE3Curve::setMinSamplingPeriod(Time time) {
  slerpPolicy_.setMinSamplingPeriod(time);
}

///   eg. 4 will add a coefficient every 4 extend
void SlerpSE3Curve::setSamplingRatio(const int ratio) {
  slerpPolicy_.setMinimumMeasurements(ratio);
}

void SlerpSE3Curve::clear() {
  manager_.clear();
}

void SlerpSE3Curve::addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const {

//  gtsam::noiseModel::Constrained::shared_ptr priorNoise = gtsam::noiseModel::Constrained::All(gtsam::traits<Coefficient>::dimension, 1e5);

  Eigen::Matrix<double,6,1> noise;
  noise(0) = 0.0000001;
  noise(1) = 0.0000001;
  noise(2) = 0.0000001;
  noise(3) = 0.0000001;
  noise(4) = 0.0000001;
  noise(5) = 0.0000001;

  gtsam::noiseModel::Diagonal::shared_ptr priorNoise = gtsam::noiseModel::Diagonal::
        Sigmas(noise);



  CoefficientIter rVal0, rVal1;
  manager_.getCoefficientsAt(priorTime, &rVal0, &rVal1);

  gtsam::ExpressionFactor<Coefficient> f0(priorNoise,
                                          rVal0->second.coefficient,
                                          gtsam::Expression<Coefficient>(rVal0->second.key));
  gtsam::ExpressionFactor<Coefficient> f1(priorNoise,
                                          rVal1->second.coefficient,
                                          gtsam::Expression<Coefficient>(rVal1->second.key));
  graph->push_back(f0);
  graph->push_back(f1);

}

void SlerpSE3Curve::transformCurve(const ValueType T) {
  std::vector<Time> coefTimes;
  manager_.getTimes(&coefTimes);
  for (size_t i = 0; i < coefTimes.size(); ++i) {
    // Apply a rigid transformation to every coefficient (on the left side).
    manager_.insertCoefficient(coefTimes[i],T*evaluate(coefTimes[i]));
  }
}

Time SlerpSE3Curve::getTimeAtKey(gtsam::Key key) const {
  return manager_.getCoefficientTimeByKey(key);
}

void SlerpSE3Curve::saveCurveTimesAndValues(const std::string& filename) const {
  std::vector<Time> curveTimes;
  manager_.getTimes(&curveTimes);

  Eigen::VectorXd v(7);

  std::vector<Eigen::VectorXd> curveValues;
  ValueType val;
  for (size_t i = 0; i < curveTimes.size(); ++i) {
    val = evaluate(curveTimes[i]);
    v << val.getPosition().x(), val.getPosition().y(), val.getPosition().z(),
        val.getRotation().w(), val.getRotation().x(), val.getRotation().y(), val.getRotation().z();
    curveValues.push_back(v);
  }

  CurvesTestHelpers::writeTimeVectorCSV(filename, curveTimes, curveValues);
}

} // namespace curves
