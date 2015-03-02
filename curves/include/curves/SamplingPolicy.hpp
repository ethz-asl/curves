/*
 * @file SamplingPolicy.hpp
 * @date Mar 02, 2015
 * @author Abel Gawel
 */

#ifndef SAMPLING_POLICY_HPP
#define SAMPLING_POLICY_HPP

namespace curves {

template<typename CurveConfig>
class SamplingPolicy {
  friend class Curve<CurveConfig>;

  typedef typename Curve<CurveConfig>::ValueType ValueType;

 protected:
  int measurementsSinceLastExtend_;
  int minimumMeasurements_;
  Time minSamplingPeriod_;

 public:

  SamplingPolicy() :
    measurementsSinceLastExtend_(0),
    minimumMeasurements_(1),
    minSamplingPeriod_(0) {};

  SamplingPolicy(int minimumMeasurements, Time minSamplingPeriod) :
    minimumMeasurements_(minimumMeasurements),
    measurementsSinceLastExtend_(0),
    minSamplingPeriod_(minSamplingPeriod) {};

  ~SamplingPolicy() {};

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
