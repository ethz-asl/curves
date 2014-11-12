#include <curves/LinearInterpolationVectorSpaceCurve.hpp>
#include <iostream>
#include "gtsam_unstable/nonlinear/Expression.h"
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace curves {

template<int N>
LinearInterpolationVectorSpaceCurve<N>::LinearInterpolationVectorSpaceCurve() : VectorSpaceCurve<N>() {}

template<int N>
LinearInterpolationVectorSpaceCurve<N>::~LinearInterpolationVectorSpaceCurve() {}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=======LINEAR INTERPOLATION CURVE========" << std::endl;
  std::cout << str << std::endl;
  std::cout << "num of coefficients: " << manager_.size() << std::endl;
  std::cout << "dimension: " << N << std::endl;
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
    manager_.getCoefficientByKey(keys[i]).print(ss.str());
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::getCoefficientsAt(const Time& time,
                                                               Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  (*outCoefficients)[rval0->key] = rval0->coefficient;
  (*outCoefficients)[rval1->key] = rval1->coefficient;
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::getCoefficientsAt(const Time& time,
                                                               KeyCoefficientTime** outCoefficient0,
                                                               KeyCoefficientTime** outCoefficient1) const {
  CHECK_NOTNULL(&outCoefficient0);
  CHECK_NOTNULL(&outCoefficient1);
  bool success = manager_.getCoefficientsAt(time, outCoefficient0, outCoefficient1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::getCoefficientsInRange(Time startTime,
                                                                    Time endTime,
                                                                    Coefficient::Map* outCoefficients) const {
  manager_.getCoefficientsInRange(startTime, endTime, outCoefficients);
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::getCoefficients(Coefficient::Map* outCoefficients) const {
  manager_.getCoefficients(outCoefficients);
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::setCoefficient(Key key, const Coefficient& value) {
  manager_.setCoefficientByKey(key, value);
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::setCoefficients(const Coefficient::Map& coefficients) {
  manager_.setCoefficients(coefficients);
}

template<int N>
Time LinearInterpolationVectorSpaceCurve<N>::getMaxTime() const {
  return manager_.getMaxTime();
}

template<int N>
Time LinearInterpolationVectorSpaceCurve<N>::getMinTime() const {
  return manager_.getMinTime();
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::fitCurve(const std::vector<Time>& times,
                                                      const std::vector<ValueType>& values,
                                                      std::vector<Key>* outKeys) {
  CHECK_EQ(times.size(), values.size());

  if(times.size() > 0) {
    manager_.clear();
    if (outKeys != NULL) {
      outKeys->clear();
      outKeys->reserve(times.size());
    }
    std::vector<Coefficient> coefficients;
    coefficients.reserve(times.size());
    size_t vsize = values[0].size();
    for(size_t i = 0; i < values.size(); ++i) {
      CHECK_EQ(vsize, values[i].size()) << "The vectors must be uniform length.";
      coefficients.push_back(Coefficient(values[i]));
    }
    manager_.insertCoefficients(times,coefficients,outKeys);
  }
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::extend(const std::vector<Time>& times,
                                                    const std::vector<ValueType>& values) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  std::vector<Key> outKeys;
  std::vector<Coefficient> coefficients(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    coefficients[i] = Coefficient(values[i]);
  }
  manager_.insertCoefficients(times, coefficients, &outKeys);
}

template<int N>
typename LinearInterpolationVectorSpaceCurve<N>::ValueType
LinearInterpolationVectorSpaceCurve<N>::evaluate(Time time) const {
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;  

  Time dt = rval1->time - rval0->time;
  Time t = rval1->time - time;
  // Alpha goes from zero to one.
  double alpha = double(t)/double(dt);

  return alpha * rval0->coefficient.getValue() + (1.0 - alpha) * rval1->coefficient.getValue();
}

template<int N>
typename LinearInterpolationVectorSpaceCurve<N>::DerivativeType
LinearInterpolationVectorSpaceCurve<N>::evaluateDerivative(Time time,
                                                           unsigned derivativeOrder) const {

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

// Evaluation function in functional form. To be passed to the expression
template<int N>
Eigen::Matrix<double,N,1> linearInterpolation(Eigen::Matrix<double,N,1>  v1,
                                              Eigen::Matrix<double,N,1>  v2, double alpha,
                                              boost::optional<Eigen::Matrix<double,N,N> &> H1=boost::none,
                                              boost::optional<Eigen::Matrix<double,N,N> &> H2=boost::none) {

  if (H1) { *H1 = Eigen::Matrix<double,N,N>::Identity()*(1-alpha); }
  if (H2) { *H2 = Eigen::Matrix<double,N,N>::Identity()*alpha; }

  return v1*(1-alpha) + v2*alpha;
}

template<int N>
gtsam::Expression<typename LinearInterpolationVectorSpaceCurve<N>::ValueType>
LinearInterpolationVectorSpaceCurve<N>::getEvalExpression(const Time& time) const {
  typedef typename LinearInterpolationVectorSpaceCurve<N>::ValueType ValueType;
  using namespace gtsam;
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);

  Expression<ValueType> leaf1(rval0->key);
  Expression<ValueType> leaf2(rval1->key);

  double alpha = double(time - rval0->time)/double(rval1->time - rval0->time);

  Expression<ValueType> rval(boost::bind(&linearInterpolation<N>, _1, _2, alpha, _3, _4),
                             leaf1, leaf2);

  return rval;
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

} // namespace curves
