/*
 * @file SdeGaussianProcessCoefficientManager.cpp
 * @date Dec 3, 2014
 * @author Sean Anderson
 */

#include <curves/SdeGaussianProcessCoefficientManager.hpp>
#include <boost/make_shared.hpp>

namespace curves {

SdeGaussianProcessCoefficientManager::SdeGaussianProcessCoefficientManager() {
}

SdeGaussianProcessCoefficientManager::~SdeGaussianProcessCoefficientManager() {
}

void SdeGaussianProcessCoefficientManager::print(const std::string& str) const {
  // \todo (Sean)
}

void SdeGaussianProcessCoefficientManager::removeCoefficientWithKey(Key key) {
  CHECK(false) << "Not implemented";
  // todo Sean -- Remove relative prior values related to this coeff, then remove
}

void SdeGaussianProcessCoefficientManager::removeCoefficientAtTime(Time time) {
  CHECK(false) << "Not implemented";
  // todo Sean -- Remove relative prior values related to this coeff, then remove
}

void SdeGaussianProcessCoefficientManager::insertNewCoefficient(Key key, Time time, const Coefficient& coefficient) {
  // Check if new coefficient is being inserted
  // between two existing coefficients
  if(time >= getMinTime() && time <= getMaxTime()) {
    std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it = timeToCoefficient_.upper_bound(time);
    --it;
    // \todo Sean -- Remove relative prior values between it->second and (++it)->second
  }
  CoefficientManagerBase::insertNewCoefficient(key, time, coefficient);
}

boost::shared_ptr<KeyCoefficientTime> SdeGaussianProcessCoefficientManager::makeContainer(Key key, Coefficient coefficient, Time time) {
  return boost::make_shared<KeyCoefficientTimePrior>(key, coefficient, time);
}

} // namespace curves
