/*
 * @file DiscreteSE3Curve.cpp
 * @date Oct 10, 2014
 * @author Renaud Dube, Abel Gawel
 */

#include <curves/DiscreteSE3Curve.hpp>
#include <iostream>

#include "gtsam/nonlinear/ExpressionFactor.h"

namespace curves {

DiscreteSE3Curve::DiscreteSE3Curve() : SE3Curve() {}

DiscreteSE3Curve::~DiscreteSE3Curve() {}

void DiscreteSE3Curve::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=========== Discrete SE3 CURVE =============" << std::endl;
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

Time DiscreteSE3Curve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time DiscreteSE3Curve::getMinTime() const {
  return manager_.getMinTime();
}

bool DiscreteSE3Curve::isEmpty() const {
  return manager_.empty();
}

int DiscreteSE3Curve::size() const {
  return manager_.size();
}

void DiscreteSE3Curve::fitCurve(const std::vector<Time>& times,
                             const std::vector<ValueType>& values,
                             std::vector<Key>* outKeys) {
  CHECK_EQ(times.size(), values.size());
  if(times.size() > 0) {
    clear();
    manager_.insertCoefficients(times,values, outKeys);
  }
}

void DiscreteSE3Curve::setCurve(const std::vector<Time>& times,
              const std::vector<ValueType>& values) {
  CHECK_EQ(times.size(), values.size());
  if(times.size() > 0) {
    manager_.insertCoefficients(times,values);
  }
}


void DiscreteSE3Curve::extend(const std::vector<Time>& times,
                           const std::vector<ValueType>& values,
                           std::vector<Key>* outKeys) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  discretePolicy_.extend<DiscreteSE3Curve, ValueType>(times, values, this, outKeys);
}

typename DiscreteSE3Curve::DerivativeType
DiscreteSE3Curve::evaluateDerivative(Time time,
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

/// \brief forms slerp interpolation into a binary expression with 2 leafs and binds alpha into it,
///        uses break down of expression into its operations
///        \f[ T = A(A^{-1}B)^{\alpha} \f]
gtsam::Expression<typename DiscreteSE3Curve::ValueType>
DiscreteSE3Curve::getValueExpression(const Time& time) const {
  typedef typename DiscreteSE3Curve::ValueType ValueType;
  using namespace gtsam;
  CoefficientIter a, b;
  bool success = manager_.getCoefficientsAt(time, &a, &b);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  Expression<ValueType> leaf_a(a->second.key);
  Expression<ValueType> leaf_b(b->second.key);

  // If the time is closer to a
  if (b->first - time >= time - a->first) {
    return leaf_a;
  } else {
    return leaf_b;
  }

}

gtsam::Expression<typename DiscreteSE3Curve::DerivativeType>
DiscreteSE3Curve::getDerivativeExpression(const Time& time, unsigned derivativeOrder) const {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

SE3 DiscreteSE3Curve::evaluate(Time time) const {
  // Check if the curve is only defined at this one time
  if (manager_.getMaxTime() == time && manager_.getMinTime() == time) {
    return manager_.coefficientBegin()->second.coefficient;
  } else {
    if (time == manager_.getMaxTime()) {
      // Efficient evaluation of a curve end
      return (--manager_.coefficientEnd())->second.coefficient;
    } else {
// Pure discrete
      CoefficientIter a, b;
      bool success = manager_.getCoefficientsAt(time, &a, &b);
      CHECK(success) << "Unable to get the coefficients at time " << time;
      SE3 T_W_A = a->second.coefficient;
      SE3 T_W_B = b->second.coefficient;

      // If the time is closer to a
      if (b->first - time >= time - a->first) {
        return T_W_A;
      } else {
        return T_W_B;
      }
    }
  }
}

void DiscreteSE3Curve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

/// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d DiscreteSE3Curve::evaluateAngularVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d DiscreteSE3Curve::evaluateAngularVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d DiscreteSE3Curve::evaluateLinearVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d DiscreteSE3Curve::evaluateLinearVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d DiscreteSE3Curve::evaluateTwistA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d DiscreteSE3Curve::evaluateTwistB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d DiscreteSE3Curve::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d DiscreteSE3Curve::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d DiscreteSE3Curve::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d DiscreteSE3Curve::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d DiscreteSE3Curve::evaluateDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d DiscreteSE3Curve::evaluateDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}

void DiscreteSE3Curve::initializeGTSAMValues(gtsam::KeySet keys, gtsam::Values* values) const {
  manager_.initializeGTSAMValues(keys, values);
}

void DiscreteSE3Curve::initializeGTSAMValues(gtsam::Values* values) const {
  manager_.initializeGTSAMValues(values);
}

void DiscreteSE3Curve::updateFromGTSAMValues(const gtsam::Values& values) {
  manager_.updateFromGTSAMValues(values);
}

void DiscreteSE3Curve::setMinSamplingPeriod(Time time) {
  discretePolicy_.setMinSamplingPeriod(time);
}

///   eg. 4 will add a coefficient every 4 extend
void DiscreteSE3Curve::setSamplingRatio(const int ratio) {
  discretePolicy_.setMinimumMeasurements(ratio);
}

void DiscreteSE3Curve::clear() {
  manager_.clear();
}

void DiscreteSE3Curve::addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const {
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

void DiscreteSE3Curve::transformCurve(const ValueType T) {
  std::vector<Time> coefTimes;
  manager_.getTimes(&coefTimes);
  for (size_t i = 0; i < coefTimes.size(); ++i) {
    // Apply a rigid transformation to every coefficient (on the left side).
    manager_.insertCoefficient(coefTimes[i],T*evaluate(coefTimes[i]));
  }
}

Time DiscreteSE3Curve::getTimeAtKey(gtsam::Key key) const {
  return manager_.getCoefficientTimeByKey(key);
}

void DiscreteSE3Curve::saveCurveTimesAndValues(const std::string& filename) const {
  std::vector<Time> curveTimes;
  manager_.getTimes(&curveTimes);

  saveCurveAtTimes(filename, curveTimes);
}

void DiscreteSE3Curve::saveCurveAtTimes(const std::string& filename, std::vector<Time> times) const {
  Eigen::VectorXd v(7);

  std::vector<Eigen::VectorXd> curveValues;
  ValueType val;
  for (size_t i = 0; i < times.size(); ++i) {
    val = evaluate(times[i]);
    v << val.getPosition().x(), val.getPosition().y(), val.getPosition().z(),
        val.getRotation().w(), val.getRotation().x(), val.getRotation().y(), val.getRotation().z();
    curveValues.push_back(v);
  }

  writeTimeVectorCSV(filename, times, curveValues);
}

void DiscreteSE3Curve::getCurveTimes(std::vector<Time>* outTimes) const {
  manager_.getTimes(outTimes);
}

} // namespace curves
