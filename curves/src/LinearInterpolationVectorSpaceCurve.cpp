#include <curves/LinearInterpolationVectorSpaceCurve.hpp>

namespace curves {

LinearInterpolationVectorSpaceCurve::LinearInterpolationVectorSpaceCurve(size_t dimension) :
    VectorSpaceCurve(dimension) {}

LinearInterpolationVectorSpaceCurve::~LinearInterpolationVectorSpaceCurve() {}

void LinearInterpolationVectorSpaceCurve::print(const std::string& str) const {

}
  

void LinearInterpolationVectorSpaceCurve::getCoefficientsAt(Time time, 
                                                            Coefficient::Map& outCoefficients) const {
  std::pair<KeyCoefficientTime*, KeyCoefficientTime*> rval;
  bool success = manager_.getCoefficientsAt(time, rval);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  outCoefficients[rval.first->key] = rval.first->coefficient;
  outCoefficients[rval.second->key] = rval.second->coefficient;
                                            
}

void LinearInterpolationVectorSpaceCurve::getCoefficientsInRange(Time startTime, 
                                                                 Time endTime, 
                                                                 Coefficient::Map& outCoefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::getCoefficients(Coefficient::Map& outCoefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}
  
void LinearInterpolationVectorSpaceCurve::setCoefficient(Key key, const Coefficient& value) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::setCoefficients(Coefficient::Map& coefficients) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}



Time LinearInterpolationVectorSpaceCurve::getBackTime() const {
  return manager_.getBackTime();
}
  
Time LinearInterpolationVectorSpaceCurve::getFrontTime() const {
  return manager_.getFrontTime();
}



void LinearInterpolationVectorSpaceCurve::extendFront(Time time) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::extendBack(Time time) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::retractFront(Time time) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::retractBack(Time time) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::extendFront(const std::vector<int64_t>& times,
                                                      const std::vector<Eigen::VectorXd>& values) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

void LinearInterpolationVectorSpaceCurve::extendBack(const std::vector<int64_t>& times,
                                                     const std::vector<Eigen::VectorXd>& values) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}


void LinearInterpolationVectorSpaceCurve::fitCurve(const std::vector<int64_t>& times,
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
    manager_.insertCoefficients(times,coefficients,outKeys);
  }
}



Eigen::VectorXd LinearInterpolationVectorSpaceCurve::evaluateVector(Time time) {
  std::pair<KeyCoefficientTime*, KeyCoefficientTime*> rval;
  bool success = manager_.getCoefficientsAt(time, rval);
  CHECK(success) << "Unable to get the coefficients at time " << time;  
  
  Time dt = rval.second->time - rval.first->time;
  Time t = rval.second->time - time;
  // Alpha goes from zero to one.
  double alpha = double(t)/double(dt);
  
  return alpha * rval.first->coefficient.getValue() + (1.0 - alpha) * rval.second->coefficient.getValue();
}
  
Eigen::VectorXd LinearInterpolationVectorSpaceCurve::evaluateDerivative(Time time, unsigned derivativeOrder) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";  
}

/// \brief Get an evaluator at this time
LinearInterpolationVectorSpaceCurve::EvaluatorTypePtr LinearInterpolationVectorSpaceCurve::getTypedEvaluator(Time time) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}


} // namespace curves
