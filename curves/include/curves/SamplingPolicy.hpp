/*
 * @file SamplingPolicy.hpp
 * @date Mar 02, 2015
 * @author Abel Gawel
 */

#ifndef SAMPLING_POLICY_HPP
#define SAMPLING_POLICY_HPP

#include <typeinfo>

namespace curves {

class SamplingPolicy {

 protected:
  int measurementsSinceLastExtend_;
  int minimumMeasurements_;
  Time minSamplingPeriod_;
  Time lastExtend_;

 public:

  SamplingPolicy() :
    measurementsSinceLastExtend_(0),
    minimumMeasurements_(1),
    minSamplingPeriod_(0),
    lastExtend_(0) {};

  SamplingPolicy(int minimumMeasurements, Time minSamplingPeriod) :
    measurementsSinceLastExtend_(0),
    minimumMeasurements_(minimumMeasurements),
    minSamplingPeriod_(minSamplingPeriod),
    lastExtend_(0) {};

  ~SamplingPolicy() {};

  template<typename CurveType, typename ValueType>
  Key interpolationExtend(const Time& time,
                          const ValueType& value,
                          CurveType* curve) {
    CHECK(false) << "no interpolation extend policy implemented for " <<  typeid(CurveType).name();
    return Key();
  }

  template<typename CurveType, typename ValueType>
  Key defaultExtend(const Time& time,
                    const ValueType& value,
                    CurveType* curve) {
    CHECK(false) << "no default extend policy implemented for " <<  typeid(CurveType).name();
    return Key();
  }

  template<typename CurveType, typename ValueType>
  void extend(const std::vector<Time>& times,
              const std::vector<ValueType>& values,
              CurveType* curve,
              std::vector<Key>* outKeys = NULL) {
    CHECK(false) << "no extend policy implemented for " <<  typeid(CurveType).name();
  }

  /// Print the value of the coefficient, for debugging and unit tests
  int getMeasurementsSinceLastExtend() {
    return measurementsSinceLastExtend_;
  }

  int getMinimumMeasurements() {
    return minimumMeasurements_;
  }

  Time getMinSamplingPeriod() {
    return minSamplingPeriod_;
  }

  Time getLastExtendTime() {
    return lastExtend_;
  }

  void setLastExtendTime(Time time) {
    lastExtend_ = time;
  }

  void setMinimumMeasurements(int n) {
    minimumMeasurements_ = n;
  }

  void setMinSamplingPeriod(Time minSamplingPeriod) {
    minSamplingPeriod_ = minSamplingPeriod;
  }

  void incrementMeasurementsTaken(int num) {
    measurementsSinceLastExtend_ = measurementsSinceLastExtend_ + num;
  }

  void setMeasurementsSinceLastExtend_(int num) {
    measurementsSinceLastExtend_ = num;
  }
};

} // namespace curves

#endif /* SAMPLING_POLICY_HPP */
