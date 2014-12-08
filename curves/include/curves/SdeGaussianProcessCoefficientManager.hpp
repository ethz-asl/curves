/*
 * @file SdeGaussianProcessCoefficientManager.hpp
 * @date Dec 3, 2014
 * @author Sean Anderson
 */

#ifndef CURVES_SDE_GAUSSIAN_PROCESS_COEFFICIENT_MANAGER_HPP
#define CURVES_SDE_GAUSSIAN_PROCESS_COEFFICIENT_MANAGER_HPP

#include "Support2CoefficientManager.hpp"

namespace curves {

struct KeyCoefficientTimePrior : public KeyCoefficientTime_CRTP<KeyCoefficientTimePrior> {
  /// \todo change example double to instead store
  /// the matrices and vectors required by the prior
  double extraPriorInfo;

  KeyCoefficientTimePrior(Key key, const Coefficient& coefficient, Time time, double extraPriorInfo) :
    KeyCoefficientTime_CRTP<KeyCoefficientTimePrior>(key, coefficient, time), extraPriorInfo(extraPriorInfo) {}
  KeyCoefficientTimePrior(Key key, const Coefficient& coefficient, Time time) :
    KeyCoefficientTime_CRTP<KeyCoefficientTimePrior>(key, coefficient, time) {}
  KeyCoefficientTimePrior() {}

  virtual bool equals(const KeyCoefficientTimePrior& other) const {
    return extraPriorInfo == other.extraPriorInfo && KeyCoefficientTime::equals(other);
  }
};

class SdeGaussianProcessCoefficientManager : public Support2CoefficientManager {
 public:

  /// Initialize a coefficient manager
  SdeGaussianProcessCoefficientManager();
  virtual ~SdeGaussianProcessCoefficientManager();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// \brief Remove the coefficient with this key.
  ///
  /// It is an error if the key does not exist.
  /// Relative prior information between this coefficient
  /// and others is removed.
  virtual void removeCoefficientWithKey(Key key);

  /// \brief Remove the coefficient at this time
  ///
  /// It is an error if there is no coefficient at this time.
  /// Relative prior information between this coefficient
  /// and others is removed.
  virtual void removeCoefficientAtTime(Time time);

 private:
  /// \brief insert a coefficient that does not exist yet
  virtual void insertNewCoefficient(Key key, Time time, const Coefficient& coefficient);

  /// Instantiate a new container using the structure KeyCoefficientTimePrior.
  virtual boost::shared_ptr<KeyCoefficientTime> makeContainer(Key key, Coefficient coefficient, Time time);
};

} // namespace curves

#endif /* CURVES_SDE_GAUSSIAN_PROCESS_COEFFICIENT_MANAGER_HPP */
