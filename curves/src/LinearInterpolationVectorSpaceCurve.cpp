#include <curves/LinearInterpolationVectorSpaceCurve.hpp>
#include <iostream>

namespace curves {

LinearInterpolationVectorSpaceCurve::LinearInterpolationVectorSpaceCurve(size_t dimension) :
                                            VectorSpaceCurve(dimension) {}

LinearInterpolationVectorSpaceCurve::~LinearInterpolationVectorSpaceCurve() {}

void LinearInterpolationVectorSpaceCurve::print(const std::string& str) const {
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


void LinearInterpolationVectorSpaceCurve::getCoefficientsAt(const Time& time,
                                                            Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  KeyCoefficientTime rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  (*outCoefficients)[rval0.key] = rval0.coefficient;
  (*outCoefficients)[rval1.key] = rval1.coefficient;
}

void LinearInterpolationVectorSpaceCurve::getCoefficientsAt(const Time& time,
                                                            KeyCoefficientTime* outCoefficient0, KeyCoefficientTime* outCoefficient1) const {
  CHECK_NOTNULL(&outCoefficient0);
  CHECK_NOTNULL(&outCoefficient1);
  if (time == this->getMaxTime()) {
    std::cout <<"max time reached" <<std::endl;
  }
  bool success = manager_.getCoefficientsAt(time, outCoefficient0, outCoefficient1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
}

void LinearInterpolationVectorSpaceCurve::getCoefficientsInRange(Time startTime, 
                                                                 Time endTime, 
                                                                 Coefficient::Map* outCoefficients) const {
  manager_.getCoefficientsInRange(startTime, endTime, outCoefficients);
}

void LinearInterpolationVectorSpaceCurve::getCoefficients(Coefficient::Map* outCoefficients) const {
  manager_.getCoefficients(outCoefficients);
}

void LinearInterpolationVectorSpaceCurve::setCoefficient(Key key, const Coefficient& value) {
  manager_.setCoefficientByKey(key, value);
}

void LinearInterpolationVectorSpaceCurve::setCoefficients(const Coefficient::Map& coefficients) {
  manager_.setCoefficients(coefficients);
}



Time LinearInterpolationVectorSpaceCurve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time LinearInterpolationVectorSpaceCurve::getMinTime() const {
  return manager_.getMinTime();
}




void LinearInterpolationVectorSpaceCurve::fitCurve(const std::vector<Time>& times,
                                                   const std::vector<Eigen::VectorXd>& values) {
  CHECK_EQ(times.size(), values.size());

  if(times.size() > 0) {
    manager_.clear();
    std::vector<Key> outKeys;
    outKeys.reserve(times.size());
    std::vector<Coefficient> coefficients;
    coefficients.reserve(times.size());
    size_t vsize = values[0].size();
    for(size_t i = 0; i < values.size(); ++i) {
      CHECK_EQ(vsize, values[i].size()) << "The vectors must be uniform length.";
      coefficients.push_back(Coefficient(values[i]));
    }
    manager_.insertCoefficients(times,coefficients,&outKeys);
  }
}

void LinearInterpolationVectorSpaceCurve::extend(const std::vector<Time>& times,
                                                 const std::vector<ValueType>& values) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  std::vector<Key> outKeys;
  std::vector<Coefficient> coefficients(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    coefficients[i] = Coefficient(values[i]);
  }
  manager_.insertCoefficients(times, coefficients, &outKeys);
}

Eigen::VectorXd LinearInterpolationVectorSpaceCurve::evaluate(Time time) const {
  KeyCoefficientTime rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;  

  Time dt = rval1.time - rval0.time;
  Time t = rval1.time - time;
  // Alpha goes from zero to one.
  double alpha = double(t)/double(dt);

  return alpha * rval0.coefficient.getValue() + (1.0 - alpha) * rval1.coefficient.getValue();
}

Eigen::VectorXd LinearInterpolationVectorSpaceCurve::evaluateDerivative(Time time, unsigned derivativeOrder) const {

  // time is out of bound --> error
  CHECK_GE(time, this->getMinTime()) << "Time out of bounds"; 
  CHECK_LE(time, this->getMaxTime()) << "Time out of bounds";

  Eigen::VectorXd dCoeff;
  Time dt;
  KeyCoefficientTime rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  // first derivative
  if (derivativeOrder == 1) {
    dCoeff = rval1.coefficient.getValue() - rval0.coefficient.getValue();
    dt = rval1.time - rval0.time;
    return dCoeff/dt;
  } else { // order of derivative > 1 returns vector of zeros
    const int dimension = rval0.coefficient.dim();
    return Eigen::VectorXd::Zero(dimension,1);
  }
}

/// \brief Get an evaluator at this time
LinearInterpolationVectorSpaceCurve::EvaluatorTypePtr LinearInterpolationVectorSpaceCurve::getEvaluator(const Time& time) const {
  boost::shared_ptr< Evaluator<VectorSpaceConfig> > rval( new LinearInterpolationVectorSpaceEvaluator((*this), time) );
  return rval;
}

void LinearInterpolationVectorSpaceCurve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

boost::unordered_map<Key, KeyCoefficientTime> LinearInterpolationVectorSpaceCurve::getKeyCoefficientTime() const {
  return manager_.getKeyCoefficientTime();
}

} // namespace curves
