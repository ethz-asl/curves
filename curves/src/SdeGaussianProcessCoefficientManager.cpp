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

boost::shared_ptr<KeyCoefficientTime> SdeGaussianProcessCoefficientManager::instantiateNewContainer(Key key, Coefficient coefficient, Time time) {
  return boost::make_shared<KeyCoefficientTimePrior>(key, coefficient, time);
}

} // namespace curves
