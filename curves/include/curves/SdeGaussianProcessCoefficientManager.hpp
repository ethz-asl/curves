/*
 * @file SdeGaussianProcessCoefficientManager.hpp
 * @date Dec 3, 2014
 * @author Sean Anderson
 */

#ifndef CURVES_SDE_GAUSSIAN_PROCESS_COEFFICIENT_MANAGER_HPP
#define CURVES_SDE_GAUSSIAN_PROCESS_COEFFICIENT_MANAGER_HPP

#include "Support2CoefficientManager.hpp"

namespace curves {

struct KeyCoefficientTimePrior : public KeyCoefficientTime {
  double extraPriorInfo;

  KeyCoefficientTimePrior(Key key, const Coefficient& coefficient, Time time, double extraPriorInfo) :
    KeyCoefficientTime(key, coefficient, time), extraPriorInfo(extraPriorInfo) {}
  KeyCoefficientTimePrior(Key key, const Coefficient& coefficient, Time time) :
    KeyCoefficientTime(key, coefficient, time) {}
  KeyCoefficientTimePrior() {}

  virtual bool equals(const KeyCoefficientTimePrior& other) const {
    return extraPriorInfo == other.extraPriorInfo && KeyCoefficientTime::equals(other);
  }
};

class SdeGaussianProcessCoefficientManager : public Support2CoefficientManager {
 public:

  SdeGaussianProcessCoefficientManager();
  virtual ~SdeGaussianProcessCoefficientManager();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

 private:
  /// Instantiate a new container using the structure KeyCoefficientTimePrior.
  virtual boost::shared_ptr<KeyCoefficientTime> instantiateNewContainer(Key key, Coefficient coefficient, Time time);
};

} // namespace curves

#endif /* CURVES_SDE_GAUSSIAN_PROCESS_COEFFICIENT_MANAGER_HPP */
