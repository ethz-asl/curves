#include <curves/GaussianProcessVectorSpaceCurve.hpp>
#include <iostream>

namespace curves {

GaussianProcessVectorSpaceCurve::GaussianProcessVectorSpaceCurve(boost::shared_ptr<GaussianProcessVectorSpacePrior> prior) :
  VectorSpaceCurve(prior->dim()), prior_(prior) { }

GaussianProcessVectorSpaceCurve::~GaussianProcessVectorSpaceCurve() {}

void GaussianProcessVectorSpaceCurve::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=======GAUSSIAN PROCESS CURVE========" << std::endl;
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

void GaussianProcessVectorSpaceCurve::getCoefficientsAt(const Time& time,
                                                        Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  (*outCoefficients)[rval0->key] = rval0->coefficient;
  (*outCoefficients)[rval1->key] = rval1->coefficient;
}

void GaussianProcessVectorSpaceCurve::getCoefficientsInRange(Time startTime,
                                                             Time endTime,
                                                             Coefficient::Map* outCoefficients) const {
  manager_.getCoefficientsInRange(startTime, endTime, outCoefficients);
}

void GaussianProcessVectorSpaceCurve::getCoefficients(Coefficient::Map* outCoefficients) const {
  manager_.getCoefficients(outCoefficients);
}

void GaussianProcessVectorSpaceCurve::setCoefficient(Key key, const Coefficient& value) {
  manager_.setCoefficientByKey(key, value);
}

void GaussianProcessVectorSpaceCurve::setCoefficients(const Coefficient::Map& coefficients) {
  manager_.setCoefficients(coefficients);
}

Time GaussianProcessVectorSpaceCurve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time GaussianProcessVectorSpaceCurve::getMinTime() const {
  return manager_.getMinTime();
}

void GaussianProcessVectorSpaceCurve::fitCurve(const std::vector<Time>& times,
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
    prior_->clearKeyTimes();
    prior_->addKeyTimes(times);
  }
}

void GaussianProcessVectorSpaceCurve::extend(const std::vector<Time>& times,
                                             const std::vector<ValueType>& values) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  std::vector<Key> outKeys;
  std::vector<Coefficient> coefficients(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    coefficients[i] = Coefficient(values[i]);
  }
  manager_.insertCoefficients(times, coefficients, &outKeys);
  prior_->addKeyTimes(times);
}

Eigen::VectorXd GaussianProcessVectorSpaceCurve::evaluate(Time time) const {
  // Get times and values of interpolation points
  // \todo replace this with a proper N-size grab of keytimes from general manager base pointer
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  std::vector<KeyCoefficientTime*> interpCoefficients;
  interpCoefficients.push_back(rval0);
  interpCoefficients.push_back(rval1);

  // Allocate memory for outputs
  std::vector<Time> keyTimes;
  std::vector<Eigen::VectorXd*> outMeanAtKeyTimes;
  std::vector<Eigen::MatrixXd*> outKtKinvMats;
  for (unsigned int i = 0; i < interpCoefficients.size(); i++) {
    keyTimes.push_back(interpCoefficients.at(i)->time);
    outMeanAtKeyTimes.push_back( new Eigen::VectorXd(this->dim()) );
    outKtKinvMats.push_back( new Eigen::MatrixXd(this->dim(),this->dim()) );
  }

  // Get interpolation information from prior
  Eigen::VectorXd mean = prior_->evaluateAndInterpMatrices(time, keyTimes, outMeanAtKeyTimes, outKtKinvMats);

  // Calculate the interpolation
  Eigen::VectorXd result = mean;
  for (unsigned int i = 0; i < interpCoefficients.size(); i++) {
    result = result + (*outKtKinvMats.at(i)) * (interpCoefficients.at(i)->coefficient.getValue() - *(outMeanAtKeyTimes.at(i)));
  }
  return result;
}

Eigen::VectorXd GaussianProcessVectorSpaceCurve::evaluateDerivative(Time time, unsigned derivativeOrder) const {
  // \todo Sean
  CHECK(false) << "Not implemented";
  return Eigen::VectorXd::Zero(this->dim());
}

/// \brief Get an evaluator at this time
typename GaussianProcessVectorSpaceCurve::EvaluatorTypePtr GaussianProcessVectorSpaceCurve::getEvaluator(const Time& time) const {
  boost::shared_ptr< Evaluator<VectorSpaceConfig> > rval( new GaussianProcessVectorSpaceEvaluator((*this), time) );
  return rval;
}

void GaussianProcessVectorSpaceCurve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Sean
  CHECK(false) << "Not implemented";
}

boost::unordered_map<Key, KeyCoefficientTime> GaussianProcessVectorSpaceCurve::getKeyCoefficientTime() const {
  return manager_.getKeyCoefficientTime();
}

} // namespace curves
