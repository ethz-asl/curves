#include <curves/HermiteCoefficientManager.hpp>
#include <curves/KeyGenerator.hpp>
#include <glog/logging.h>

namespace curves {

HermiteCoefficientManager::HermiteCoefficientManager()  { }
HermiteCoefficientManager::~HermiteCoefficientManager() { }

/// Compare this Coeficient with another for equality.
bool HermiteCoefficientManager::equals(const HermiteCoefficientManager& other, double tol) const {
  bool equal = true;
  equal &= coefficients_.size() == other.coefficients_.size();
  equal &= timeToCoefficient_.size() == other.timeToCoefficient_.size();
  if(equal) {
    boost::unordered_map<Key, KeyCoefficientTime>::const_iterator it1, it2;
    it1 = coefficients_.begin();
    it2 = other.coefficients_.begin();
    for( ; it1 != coefficients_.end(); ++it1, ++it2) {
      equal &= it1->first == it2->first;
      equal &= it1->second.equals(it2->second);
    }

    std::map<Time, KeyCoefficientTime*>::const_iterator it3, it4;
    it3 = timeToCoefficient_.begin();
    it4 = other.timeToCoefficient_.begin();
    for( ; it3 != timeToCoefficient_.end(); ++it3, ++it4) {
      equal &= it3->first == it4->first;
      if(it3->second && it4->second) {
        equal &= it3->second->equals(*(it4->second));
      } else {
        equal &= it3->second == it4->second;
      }
    }
    
  }
  return equal;
}
 
void HermiteCoefficientManager::getKeys(std::vector<Key>& outKeys) const {
  outKeys.clear();
  appendKeys(outKeys);
}
void HermiteCoefficientManager::appendKeys(std::vector<Key>& outKeys) const {
  outKeys.reserve(outKeys.size() + coefficients_.size());
  std::map<Time, KeyCoefficientTime*>::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    outKeys.push_back(it->second->key);
  }
}

void HermiteCoefficientManager::getTimes(std::vector<Time>& outTimes) const {
  outTimes.clear();
  outTimes.reserve(timeToCoefficient_.size());
  std::map<Time, KeyCoefficientTime*>::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    outTimes.push_back(it->first);
  }
}

void HermiteCoefficientManager::print(const std::string& str) const {
  // \todo (Abel or Renaud)
}

Key HermiteCoefficientManager::insertCoefficient(Time time, const Coefficient& coefficient) {
  typedef boost::unordered_map<Key, KeyCoefficientTime>::iterator iterator;
  Key key = KeyGenerator::getNextKey();
  std::pair<iterator, bool > success = 
      coefficients_.insert( 
          std::make_pair(key, KeyCoefficientTime(key, coefficient, time)));
  timeToCoefficient_[time] = &(success.first->second);
  return key;
}

/// \brief insert coefficients. Returns the keys for these coefficients
void HermiteCoefficientManager::insertCoefficients(const std::vector<Time>& times,
                        const std::vector<Coefficient>& values,
                        std::vector<Key>& outKeys) {
  CHECK_EQ(times.size(), values.size());
  for(Key i = 0; i < times.size(); ++i) {
    outKeys.push_back(insertCoefficient(times[i], values[i]));
  }
}
  
/// \brief return true if there is a coefficient at this time
bool HermiteCoefficientManager::hasCoefficientAtTime(Time time) const {
  std::map<Time, KeyCoefficientTime*>::const_iterator it = timeToCoefficient_.find(time);
  return it != timeToCoefficient_.end();
}

/// \brief return true if there is a coefficient with this key
bool HermiteCoefficientManager::hasCoefficientWithKey(Key key) const {
  boost::unordered_map<Key, KeyCoefficientTime>::const_iterator it = coefficients_.find(key);
  return it != coefficients_.end();
}

/// \brief set the coefficient associated with this key
///
/// This function fails if there is no coefficient associated
/// with this key.
void HermiteCoefficientManager::setCoefficientByKey(Key key, const Coefficient& coefficient) {
  boost::unordered_map<Key, KeyCoefficientTime>::iterator it = coefficients_.find(key);
  CHECK( it != coefficients_.end() ) << "Key " << key << " is not in the container.";
  it->second.coefficient = coefficient;
  
}

/// \brief set the coefficient associated with this key
///
/// This function fails if there is no coefficient associated
/// with this key.
void HermiteCoefficientManager::setCoefficientVectorByKey(Key key, const Eigen::VectorXd& vector) {
  boost::unordered_map<Key, KeyCoefficientTime>::iterator it = coefficients_.find(key);
  CHECK( it != coefficients_.end() ) << "Key " << key << " is not in the container.";
  it->second.coefficient.setVector(vector);
}

/// \brief get the coefficient associated with this key
Coefficient HermiteCoefficientManager::getCoefficientByKey(Key key) const {
  boost::unordered_map<Key, KeyCoefficientTime>::const_iterator it = coefficients_.find(key);
  CHECK( it != coefficients_.end() ) << "Key " << key << " is not in the container.";
  return it->second.coefficient;
}
  

/// \brief Get the coefficients that are active at a certain time.
bool HermiteCoefficientManager::getCoefficientsAt(Time time, 
  std::pair<KeyCoefficientTime*, KeyCoefficientTime*>& outCoefficients) const {
  
  if( timeToCoefficient_.empty() ) {
    LOG(INFO) << "No coefficients";
    return false;
  }

  std::map<Time, KeyCoefficientTime*>::const_iterator it;
  it = timeToCoefficient_.upper_bound(time);
  if(time == getFrontTime()) {
    --it;
  }
  if(it == timeToCoefficient_.begin() || it == timeToCoefficient_.end()) {
    LOG(INFO) << "time, " << time << ", is out of bounds: [" << getBackTime() << ", " << getFrontTime() << "]";
    return false;
  }
  --it;

  // Okay. Now we know that the time is bracketed by
  // it and it + 1.
  outCoefficients.first = it->second;
  outCoefficients.second = (++it)->second;
  
  return true;
}
  
/// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
void HermiteCoefficientManager::getCoefficientsInRange(Time startTime, 
                            Time endTime, 
                            Coefficient::Map& outCoefficients) const {
  std::map<Time, KeyCoefficientTime*>::const_iterator it;
  it = timeToCoefficient_.lower_bound(startTime);
  for( ; it != timeToCoefficient_.end() && it->first < endTime; ++it) {
    outCoefficients[it->second->key] = it->second->coefficient;
  }
  if( it != timeToCoefficient_.end() ) {
    outCoefficients[it->second->key] = it->second->coefficient;
  }
  
}
  
/// \brief Get all of the curve's coefficients.
void HermiteCoefficientManager::getCoefficients(Coefficient::Map& outCoefficients) const {
  std::map<Time, KeyCoefficientTime*>::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    outCoefficients[it->second->key] = it->second->coefficient;
  }
}

/// \brief Set coefficients.
///
/// If any of these coefficients doen't exist, there is an error
void HermiteCoefficientManager::setCoefficients(Coefficient::Map& coefficients) {
  // \todo Abel or Renaud
  CHECK(false) << "Not implemented";
}

/// \brief return the number of coefficients
Key HermiteCoefficientManager::size() const {
  return coefficients_.size();
}
  
/// \brief clear the coefficients
void HermiteCoefficientManager::clear() {
  coefficients_.clear();
  timeToCoefficient_.clear();
}


Time HermiteCoefficientManager::getBackTime() const {
  if(timeToCoefficient_.empty()) {
    return 0;
  }
  return timeToCoefficient_.begin()->first;
}

Time HermiteCoefficientManager::getFrontTime() const {
  if(timeToCoefficient_.empty()) {
    return 0;
  }
  return timeToCoefficient_.rbegin()->first;
}

void HermiteCoefficientManager::checkInternalConsistency(bool doExit) const {
  CHECK_EQ(coefficients_.size(), timeToCoefficient_.size());
  std::map<Time, KeyCoefficientTime*>::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    CHECK_NOTNULL(it->second);
    CHECK_EQ(it->first, it->second->time);
    Key key = it->second->key;
    boost::unordered_map<Key, KeyCoefficientTime>::const_iterator itc = coefficients_.find(key);
    CHECK( itc != coefficients_.end() ) << "Key " << key << " is not in the map";
    // This is probably the important one.
    // Check that the it->second pointer
    // points to the same object as itc->second.
    // It is supposedly guaranteed by the 
    // unordered map interface that these
    // pointers never get reallocated.
    CHECK_EQ( &itc->second, it->second);
  }
  if( doExit ) {
    exit(0);
  }
}

void HermiteCoefficientManager::removeCoefficientWithKey(Key key) {
  CHECK(false) << "Not implemented";
  // todo Abel and Renaud
}

void HermiteCoefficientManager::removeCoefficientAtTime(Time time) {
  CHECK(false) << "Not implemented";
  // todo Abel and Renaud
}


} // namespace curves


