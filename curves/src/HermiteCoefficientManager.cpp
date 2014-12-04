/*
 * @file HermiteCoefficientManager.cpp
 * @date Aug 17, 2014
 * @author Paul Furgale, Abel Gawel, Renaud Dube
 */

#include <curves/HermiteCoefficientManager.hpp>
#include <curves/KeyGenerator.hpp>
#include <glog/logging.h>

namespace curves {

HermiteCoefficientManager::HermiteCoefficientManager() {
}
HermiteCoefficientManager::~HermiteCoefficientManager() {
}

void HermiteCoefficientManager::print(const std::string& str) const {
  // \todo (Abel or Renaud)
}

/// \brief Get the coefficients that are active at a certain time.
bool HermiteCoefficientManager::getCoefficientsAt(Time time, 
                                                  KeyCoefficientTime** outCoefficient0,
                                                  KeyCoefficientTime** outCoefficient1) const {
  CHECK_NOTNULL(outCoefficient0);
  CHECK_NOTNULL(outCoefficient1);
  if( timeToCoefficient_.empty() ) {
    LOG(INFO) << "No coefficients";
    return false;
  }

  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it;

  if(time == getMaxTime()) {
    it = timeToCoefficient_.end();
    --it;
  } else {
    it = timeToCoefficient_.upper_bound(time);
  }
  if(it == timeToCoefficient_.begin() || it == timeToCoefficient_.end()) {
    LOG(INFO) << "time, " << time << ", is out of bounds: [" << getMinTime() << ", " << getMaxTime() << "]";
    return false;
  }
  --it;

  // Okay. Now we know that the time is bracketed by
  // it and it + 1.
  *outCoefficient0 = it->second.get();
  *outCoefficient1 = (++it)->second.get();

  return true;
}

} // namespace curves
