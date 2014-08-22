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


void LinearInterpolationVectorSpaceCurve::getCoefficientsAt(Time time, 
                                                            Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  std::pair<KeyCoefficientTime*, KeyCoefficientTime*> rval;
  bool success = manager_.getCoefficientsAt(time, &rval);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  (*outCoefficients)[rval.first->key] = rval.first->coefficient;
  (*outCoefficients)[rval.second->key] = rval.second->coefficient;
                                            
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
  std::pair<KeyCoefficientTime*, KeyCoefficientTime*> rval;
  bool success = manager_.getCoefficientsAt(time, &rval);
  CHECK(success) << "Unable to get the coefficients at time " << time;  

  Time dt = rval.second->time - rval.first->time;
  Time t = rval.second->time - time;
  // Alpha goes from zero to one.
  double alpha = double(t)/double(dt);

  return alpha * rval.first->coefficient.getValue() + (1.0 - alpha) * rval.second->coefficient.getValue();
}

Eigen::VectorXd LinearInterpolationVectorSpaceCurve::evaluateDerivative(Time time, unsigned derivativeOrder) const {

  // time is out of bound --> error
  CHECK_GE(time, this->getMinTime()) << "Time out of bounds"; 
  CHECK_LE(time, this->getMaxTime()) << "Time out of bounds";

  Eigen::VectorXd dCoeff;
  Time dt;
  std::pair<KeyCoefficientTime*, KeyCoefficientTime*> rval;
  bool success = manager_.getCoefficientsAt(time, &rval);
  // first derivative
  if (derivativeOrder == 1) {
    dCoeff = rval.second->coefficient.getValue() - rval.first->coefficient.getValue();
    dt = rval.second->time - rval.first->time;
    return dCoeff/dt;
  } else { // order of derivative > 1 returns vector of zeros
    const int dimension = rval.first->coefficient.dim();
    return Eigen::VectorXd::Zero(dimension,1);
  }
}

/// \brief Get an evaluator at this time
LinearInterpolationVectorSpaceCurve::EvaluatorTypePtr LinearInterpolationVectorSpaceCurve::getEvaluator(Time time) const {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}


} // namespace curves
