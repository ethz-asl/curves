/*
 * @file CompositionCurve-inl.hpp
 * @date Feb 06, 2015
 * @author Abel Gawel, Renaud Dubé, Mike Bosse
 */

#include "SE3CompositionCurve.hpp"

namespace curves{

template <class C1, class C2>
SE3CompositionCurve<C1, C2>::SE3CompositionCurve() {

}

template <class C1, class C2>
SE3CompositionCurve<C1, C2>::~SE3CompositionCurve() {

}
template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::print(const std::string& str) const {
  baseCurve_.print("Base curve");
  correctionCurve_.print("Correction curve");
}

template <class C1, class C2>
Time SE3CompositionCurve<C1, C2>::getMinTime() const{
  return baseCurve_.getMinTime();
}

template <class C1, class C2>
Time SE3CompositionCurve<C1, C2>::getMaxTime() const{
  return baseCurve_.getMaxTime();
}

template <class C1, class C2>
bool SE3CompositionCurve<C1, C2>::isEmpty() const{
  return baseCurve_.isEmpty();
}

template <class C1, class C2>
int SE3CompositionCurve<C1, C2>::size() const{
  return baseCurve_.size();
}

template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::extend(const std::vector<Time>& times,
                                         const std::vector<typename SE3CompositionCurve<C1, C2>::ValueType>& values,
                                         std::vector<Key>* outKeys){
  baseCurve_.extend(times, values, outKeys);

  // Make sure that the correction curve can be evaluated at all times
  // where the base curve is defined.
  std::vector<Time> correctionTimes;
  std::vector<ValueType> correctionValues;
  if(correctionCurve_.getMinTime() > baseCurve_.getMinTime()) {
    correctionTimes.push_back(baseCurve_.getMinTime());
    ValueType identityTransform;
    identityTransform.setIdentity();
    correctionValues.push_back(identityTransform);
  }
  if(correctionCurve_.getMaxTime() < baseCurve_.getMaxTime()) {
    correctionTimes.push_back(baseCurve_.getMaxTime());
    ValueType identityTransform;
    identityTransform.setIdentity();
    correctionValues.push_back(identityTransform);
  }
  correctionCurve_.extend(correctionTimes, correctionValues);
}

template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::fitCurve(const std::vector<Time>& times,
                                           const std::vector<typename SE3CompositionCurve<C1, C2>::ValueType>& values,
                                           std::vector<Key>* outKeys){
  extend(times, values, outKeys);
}
template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::addCorrectionCoefficients(const std::vector<Time>& times) {
  std::vector<ValueType> values;
  ValueType identityTransform;
  identityTransform.setIdentity();
  for (size_t i = 0; i < times.size(); ++i) {
    values.push_back(identityTransform);
  }
  correctionCurve_.extend(times, values);
}

template <class C1, class C2>
typename SE3CompositionCurve<C1, C2>::ValueType SE3CompositionCurve<C1, C2>::evaluate(Time time) const{
  return correctionCurve_.evaluate(time) * baseCurve_.evaluate(time);
}

template <class C1, class C2>
typename SE3CompositionCurve<C1, C2>::DerivativeType SE3CompositionCurve<C1, C2>::evaluateDerivative(Time time,
                                                                                                     unsigned derivativeOrder) const{
  //todo
}

template <class C1, class C2>
gtsam::Expression<typename SE3CompositionCurve<C1, C2>::ValueType>
SE3CompositionCurve<C1, C2>::getValueExpression(const Time& time) const{
  return kindr::minimal::compose(correctionCurve_.getValueExpression(time),
                                 gtsam::Expression<ValueType>(baseCurve_.evaluate(time)));
}

template <class C1, class C2>
gtsam::Expression<typename SE3CompositionCurve<C1, C2>::DerivativeType>
SE3CompositionCurve<C1, C2>::getDerivativeExpression(const Time& time, unsigned derivativeOrder) const{
  //todo
}

template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::setTimeRange(Time minTime, Time maxTime){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateAngularVelocityA(Time time){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateAngularVelocityB(Time time){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateLinearVelocityA(Time time){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateLinearVelocityB(Time time){
  //todo
}

template <class C1, class C2>
Vector6d SE3CompositionCurve<C1, C2>::evaluateTwistA(Time time){
  //todo
}

template <class C1, class C2>
Vector6d SE3CompositionCurve<C1, C2>::evaluateTwistB(Time time){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time){
  //todo
}

template <class C1, class C2>
Eigen::Vector3d SE3CompositionCurve<C1, C2>::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time){
  //todo
}

template <class C1, class C2>
Vector6d SE3CompositionCurve<C1, C2>::evaluateDerivativeA(unsigned derivativeOrder, Time time){
  //todo
}

template <class C1, class C2>
Vector6d SE3CompositionCurve<C1, C2>::evaluateDerivativeB(unsigned derivativeOrder, Time time){
  //todo
}

template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) const{
  correctionCurve_.initializeGTSAMValues(keys, values);
}

template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::initializeGTSAMValues(gtsam::Values* values) const{
  correctionCurve_.initializeGTSAMValues(values);
}

template <class C1, class C2>
void SE3CompositionCurve<C1, C2>::updateFromGTSAMValues(const gtsam::Values& values){
  correctionCurve_.updateFromGTSAMValues(values);
}

} // namespace curves
