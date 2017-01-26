/*

 * @file LinearInterpolationVectorSpaceCurve.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale, Abel Gawel, Renaud Dube


#include <curves/LinearInterpolationVectorSpaceCurve.hpp>
#include <iostream>
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
  std::cout << "number of coefficients: " << manager_.size() << std::endl;
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
    std::cout << ss.str() << manager_.getCoefficientByKey(keys[i]).transpose() << std::endl;
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
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
                                                    const std::vector<ValueType>& values,
                                                    std::vector<Key>* outKeys) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  std::vector<Coefficient> coefficients(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    coefficients[i] = Coefficient(values[i]);
  }
  manager_.insertCoefficients(times, coefficients, outKeys);
}

template<int N>
typename LinearInterpolationVectorSpaceCurve<N>::ValueType
LinearInterpolationVectorSpaceCurve<N>::evaluate(Time time) const {
  CoefficientIter rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;

  Time dt = rval1->first - rval0->first;
  Time t = rval1->first - time;
  // Alpha goes from zero to one.
  double alpha = double(t)/double(dt);

  return alpha * rval0->second.coefficient + (1.0 - alpha) * rval1->second.coefficient;
}

template<int N>
typename LinearInterpolationVectorSpaceCurve<N>::DerivativeType
LinearInterpolationVectorSpaceCurve<N>::evaluateDerivative(Time time,
                                                           unsigned derivativeOrder) const {

  // time is out of bound --> error
  CHECK_GE(time, this->getMinTime()) << "Time out of bounds";
  CHECK_LT(time, this->getMaxTime()) << "Time out of bounds";
  CHECK_GT(derivativeOrder, 0) << "DerivativeOrder must be greater than 0";

  typename LinearInterpolationVectorSpaceCurve<N>::DerivativeType dCoeff;
  Time dt;
  CoefficientIter rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  // first derivative
  if (derivativeOrder == 1) {
    dCoeff = rval1->second.coefficient - rval0->second.coefficient;
    dt = rval1->first - rval0->first;
    return dCoeff/dt;
  } else { // order of derivative > 1 returns vector of zeros
    dCoeff.Zero();
    return dCoeff;
  }
}

// Evaluation function in functional form. To be passed to the expression
template<int N>
Eigen::Matrix<double,N,1> linearInterpolation(Eigen::Matrix<double,N,1>  v1,
                                              Eigen::Matrix<double,N,1>  v2, double alpha,
                                              gtsam::OptionalJacobian<N,N> H1,
                                              gtsam::OptionalJacobian<N,N> H2) {
  if (H1) { *H1 = Eigen::Matrix<double,N,N>::Identity()*(1-alpha); }
  if (H2) { *H2 = Eigen::Matrix<double,N,N>::Identity()*alpha; }

  return v1*(1-alpha) + v2*alpha;
}

template<int N>
gtsam::Expression<typename LinearInterpolationVectorSpaceCurve<N>::ValueType>
LinearInterpolationVectorSpaceCurve<N>::getValueExpression(const Time& time) const {
  typedef typename LinearInterpolationVectorSpaceCurve<N>::ValueType ValueType;
  using namespace gtsam;
  CoefficientIter rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;

  Expression<ValueType> leaf1(rval0->second.key);
  Expression<ValueType> leaf2(rval1->second.key);

  double alpha = double(time - rval0->first)/double(rval1->first - rval0->first);

  Expression<ValueType> rval(boost::bind(&linearInterpolation<N>, _1, _2, alpha, _3, _4),
                             leaf1, leaf2);

  return rval;
}

template<int N>
gtsam::Expression<typename LinearInterpolationVectorSpaceCurve<N>::DerivativeType>
LinearInterpolationVectorSpaceCurve<N>::getDerivativeExpression(const Time& time, unsigned derivativeOrder) const {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::initializeGTSAMValues(gtsam::KeySet keys, gtsam::Values* values) const {
  manager_.initializeGTSAMValues(keys, values);
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::initializeGTSAMValues(gtsam::Values* values) const {
  manager_.initializeGTSAMValues(values);
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::updateFromGTSAMValues(const gtsam::Values& values) {
  manager_.updateFromGTSAMValues(values);
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::clear() {
  manager_.clear();
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const {
  gtsam::noiseModel::Constrained::shared_ptr priorNoise = gtsam::noiseModel::Constrained::All(gtsam::traits<Coefficient>::dimension);

  // Constraint the coefficients which influence the curve value at priorTime
  CoefficientIter rVal0, rVal1;
  manager_.getCoefficientsAt(priorTime, &rVal0, &rVal1);

  gtsam::ExpressionFactor<Coefficient> factor0(priorNoise,
                                          rVal0->second.coefficient,
                                          gtsam::Expression<Coefficient>(rVal0->second.key));
  gtsam::ExpressionFactor<Coefficient> factor1(priorNoise,
                                          rVal1->second.coefficient,
                                          gtsam::Expression<Coefficient>(rVal1->second.key));
  graph->push_back(factor0);
  graph->push_back(factor1);
}

template<int N>
void LinearInterpolationVectorSpaceCurve<N>::transformCurve(const ValueType T) {
  //todo
}

template<int N>
Time LinearInterpolationVectorSpaceCurve<N>::getTimeAtKey(gtsam::Key key) const {
  return manager_.getCoefficientTimeByKey(key);
}

} // namespace curves
*/
