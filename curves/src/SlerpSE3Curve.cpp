#include <curves/SlerpSE3Curve.hpp>
#include <curves/SE3CoefficientImplementation.hpp>
#include <iostream>

namespace curves {

SlerpSE3Curve::SlerpSE3Curve() : SE3Curve() {}

SlerpSE3Curve::~SlerpSE3Curve() {}

void SlerpSE3Curve::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=======LINEAR INTERPOLATION CURVE========" << std::endl;
  std::cout << str << std::endl;
  std::cout << "num of coefficients: " << manager_.size() << std::endl;
  std::cout << "dimension: " << dim() << std::endl;
  std::stringstream ss;
  std::vector<Key> keys;
  std::vector<Time> times;
  manager_.getTimes(&times);
  manager_.getKeys(&keys);
  std::cout << "curve defined between times: " << manager_.getMinTime() << " and " << manager_.getMaxTime() <<std::endl;
  std::cout <<"=========================================" <<std::endl;
  for (size_t i = 0; i < manager_.size(); i++) {
    ss << "coefficient " << keys[i] << ": ";
    manager_.getCoefficientByKey(keys[i]).print(ss.str());
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
}


void SlerpSE3Curve::getCoefficientsAt(const Time& time,
                                                            Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  (*outCoefficients)[rval0->key] = rval0->coefficient;
  (*outCoefficients)[rval1->key] = rval1->coefficient;
}

void SlerpSE3Curve::getCoefficientsAt(const Time& time,
                                                            KeyCoefficientTime** outCoefficient0,
                                                            KeyCoefficientTime** outCoefficient1) const {
  CHECK_NOTNULL(&outCoefficient0);
  CHECK_NOTNULL(&outCoefficient1);
  if (time == this->getMaxTime()) {
    std::cout <<"max time reached" <<std::endl;
  }
  bool success = manager_.getCoefficientsAt(time, outCoefficient0, outCoefficient1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
}

void SlerpSE3Curve::getCoefficientsInRange(Time startTime,
                                                                 Time endTime, 
                                                                 Coefficient::Map* outCoefficients) const {
  manager_.getCoefficientsInRange(startTime, endTime, outCoefficients);
}

void SlerpSE3Curve::getCoefficients(Coefficient::Map* outCoefficients) const {
  manager_.getCoefficients(outCoefficients);
}

void SlerpSE3Curve::setCoefficient(Key key, const Coefficient& value) {
  manager_.setCoefficientByKey(key, value);
}

void SlerpSE3Curve::setCoefficients(const Coefficient::Map& coefficients) {
  manager_.setCoefficients(coefficients);
}

Time SlerpSE3Curve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time SlerpSE3Curve::getMinTime() const {
  return manager_.getMinTime();
}

void SlerpSE3Curve::fitCurve(const std::vector<Time>& times,
                                                   const std::vector<SlerpSE3Curve::ValueType>& values) {
  CHECK_EQ(times.size(), values.size());

  if(times.size() > 0) {
    Eigen::VectorXd val(7);
    manager_.clear();
    std::vector<Key> outKeys;
    outKeys.reserve(times.size());
    std::vector<Coefficient> coefficients;
    coefficients.reserve(times.size());
    for(size_t i = 0; i < values.size(); ++i) {
      CoefficientImplementation::Ptr impl(new SE3CoefficientImplementation);
      boost::dynamic_pointer_cast<SE3CoefficientImplementation>(impl)->makeValue(values[i],&val);
      Coefficient c1(impl,val);
      coefficients.push_back(c1);
    }
    manager_.insertCoefficients(times,coefficients,&outKeys);
  }
}

void SlerpSE3Curve::extend(const std::vector<Time>& times,
                                                 const std::vector<ValueType>& values) {

//  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
//  std::vector<Key> outKeys;
//  std::vector<Coefficient> coefficients(values.size());
//  for (size_t i = 0; i < values.size(); ++i) {
//    coefficients[i] = Coefficient(values[i]);
//  }
//  manager_.insertCoefficients(times, coefficients, &outKeys);
}

SlerpSE3Curve::ValueType SlerpSE3Curve::evaluate(Time time) const {
//  KeyCoefficientTime *rval0, *rval1;
//  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
//  CHECK(success) << "Unable to get the coefficients at time " << time;
//
//  Time dt = rval1->time - rval0->time;
//  Time t = rval1->time - time;
//  // Alpha goes from zero to one.
//  double alpha = double(t)/double(dt);
//
//  return alpha * rval0->coefficient.getValue() + (1.0 - alpha) * rval1->coefficient.getValue();
}

SlerpSE3Curve::DerivativeType SlerpSE3Curve::evaluateDerivative(Time time, unsigned derivativeOrder) const {

  // time is out of bound --> error
  CHECK_GE(time, this->getMinTime()) << "Time out of bounds";
  CHECK_LE(time, this->getMaxTime()) << "Time out of bounds";

  Eigen::VectorXd dCoeff;
  Time dt;
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  // first derivative
  if (derivativeOrder == 1) {
    dCoeff = rval1->coefficient.getValue() - rval0->coefficient.getValue();
    dt = rval1->time - rval0->time;
    return dCoeff/dt;
  } else { // order of derivative > 1 returns vector of zeros
    const int dimension = rval0->coefficient.dim();
    return Eigen::VectorXd::Zero(dimension,1);
  }
}

/// \brief Get an evaluator at this time
SlerpSE3Curve::EvaluatorTypePtr SlerpSE3Curve::getEvaluator(const Time& time) const {
  boost::shared_ptr< Evaluator<SE3Config> > rval( new SlerpSE3Evaluator((*this), time) );
  return rval;
}

gtsam::Expression<SlerpSE3Curve::ValueType> SlerpSE3Curve::getEvalExpression(const Time& time) const {
  return gtsam::Expression<SlerpSE3Curve::ValueType>(1);
}

void SlerpSE3Curve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

/// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityA(Time time) {

}

/// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityB(Time time) {

}

/// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityA(Time time) {

}

/// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityB(Time time) {

}

/// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
///        expressed in Frame a. The return value has the linear velocity (0,1,2),
///        and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistA(Time time) {

}

/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
///        expressed in Frame b. The return value has the linear velocity (0,1,2),
///        and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistB(Time time) {

}

/// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) {

}

/// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) {

}

/// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) {

}

/// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) {

}

/// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
///        expressed in Frame a. The return value has the linear velocity (0,1,2),
///        and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeA(unsigned derivativeOrder, Time time) {

}

/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
///        expressed in Frame b. The return value has the linear velocity (0,1,2),
///        and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeB(unsigned derivativeOrder, Time time) {

}



} // namespace curves
