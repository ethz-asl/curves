/*
 * LocalSupport2CoefficientManager-inl.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Paul Furgale, Abel Gawel, Renaud Dube, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <curves/LocalSupport2CoefficientManager.hpp>

#include <iostream>
#include <curves/LocalSupport2CoefficientManager.hpp>
#include <curves/KeyGenerator.hpp>
#include <glog/logging.h>

namespace curves {

template <class Coefficient>
LocalSupport2CoefficientManager<Coefficient>::LocalSupport2CoefficientManager() {
}

template <class Coefficient>
LocalSupport2CoefficientManager<Coefficient>::~LocalSupport2CoefficientManager() {
}

/// Compare this Coefficient manager with another for equality.
template <class Coefficient>
bool LocalSupport2CoefficientManager<Coefficient>::equals(const LocalSupport2CoefficientManager& other,
                                                          double tol) const {
  bool equal = true;
  equal &= keyToCoefficient_.size() == other.keyToCoefficient_.size();
  equal &= timeToCoefficient_.size() == other.timeToCoefficient_.size();
  if (equal) {
    CoefficientIter it1, it2;
    it1 = keyToCoefficient_.begin();
    it2 = other.keyToCoefficient_.begin();
    for( ; it1 != keyToCoefficient_.end(); ++it1, ++it2) {
      equal &= it1->first == it2->first;
      equal &= it1->second.equals(it2->second);
    }

    CoefficientIter it3, it4;
    it3 = timeToCoefficient_.begin();
    it4 = other.timeToCoefficient_.begin();
    for( ; it3 != timeToCoefficient_.end(); ++it3, ++it4) {
      equal &= it3->first == it4->first;
      if (it3->second && it4->second) {
        equal &= it3->second->equals(*(it4->second));
      } else {
        equal &= it3->second == it4->second;
      }
    }

  }
  return equal;
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::getKeys(std::vector<Key>* outKeys) const {
  CHECK_NOTNULL(outKeys);
  outKeys->clear();
  appendKeys(outKeys);
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::appendKeys(std::vector<Key>* outKeys) const {
  CHECK_NOTNULL(outKeys);
  outKeys->reserve(outKeys->size() + keyToCoefficient_.size());
  CoefficientIter it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    outKeys->push_back(it->second.key);
  }
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::getTimes(std::vector<Time>* outTimes) const {
  CHECK_NOTNULL(outTimes);
  outTimes->clear();
  outTimes->reserve(timeToCoefficient_.size());
  CoefficientIter it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    outTimes->push_back(it->first);
  }
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::getTimesInWindow(std::vector<Time>* outTimes,
                                                                    Time begTime, Time endTime) const {
  CHECK_EQ(endTime, getMaxTime()) << "Not implemented for window not at the end.";
  CHECK(begTime >= getMinTime()) << "Asked for times outside the curve.";
  CHECK_NOTNULL(outTimes);

  outTimes->clear();
  CoefficientIter it = --(timeToCoefficient_.end());

  do {
    if (it->first >= begTime) {
      outTimes->push_back(it->first);
    }
    --it;
  } while (it->first >= begTime && it != timeToCoefficient_.begin());

  std::reverse(outTimes->begin(),outTimes->end());
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::print(const std::string& str) const {
  // \todo (Abel or Renaud)
}

template <class Coefficient>
Key LocalSupport2CoefficientManager<Coefficient>::insertCoefficient(Time time, const Coefficient& coefficient) {
  CoefficientIter it;
  Key key;

  if (this->hasCoefficientAtTime(time, &it)) {
    this->updateCoefficientByKey(it->second.key, coefficient);
    key = it->second.key;
  } else {
    key = KeyGenerator::getNextKey();
    std::pair<Time, KeyCoefficient> iterator(time, KeyCoefficient(key, coefficient));
    std::pair<CoefficientIter, bool> success =
        timeToCoefficient_.insert(iterator);
    keyToCoefficient_[key] = success.first;
  }
  return key;
}

/// \brief insert coefficients. Optionally returns the keys for these coefficients
template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::insertCoefficients(const std::vector<Time>& times,
                                                                      const std::vector<Coefficient>& values,
                                                                      std::vector<Key>* outKeys) {
  CHECK_EQ(times.size(), values.size());
  for(Key i = 0; i < times.size(); ++i) {
    if (outKeys != NULL) {
      outKeys->push_back(insertCoefficient(times[i], values[i]));
    } else {
      insertCoefficient(times[i], values[i]);
    }
  }
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::modifyCoefficientsValuesInBatch(const std::vector<Time>& times,
                                                                                   const std::vector<Coefficient>& values) {
  CHECK_EQ(times.size(), values.size());
  // Get an iterator to the first coefficient
  typename TimeToKeyCoefficientMap::iterator it = timeToCoefficient_.end();

  do {
    --it;
  } while (it->first != times[0]);

  for (size_t i = 0; i < times.size(); ++i) {
    CHECK_EQ(it->first,times[i]);
    it->second.coefficient = values[i];
    ++it;
  }
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::addCoefficientAtEnd(Time time, const Coefficient& coefficient, std::vector<Key>* outKeys) {
  CHECK(time > getMaxTime()) << "Time to add is not greater than curve max time";

  Key key = KeyGenerator::getNextKey();

  // Insert the coefficient with a hint that it goes at the end
  CoefficientIter it = timeToCoefficient_.insert(--(timeToCoefficient_.end()),
                                                 std::pair<Time, KeyCoefficient>(time, KeyCoefficient(key, coefficient)));

  keyToCoefficient_.insert(keyToCoefficient_.end(), std::pair<Key, CoefficientIter>(key,it));
  if (outKeys != NULL) {
    outKeys->push_back(key);
  }
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::modifyCoefficient(typename TimeToKeyCoefficientMap::iterator it,
                                                                     Time time, const Coefficient& coefficient) {
  // This is used by slerp sampling policy.
  // In this case a new coefficient should be placed slightly later than the initial one.
  CoefficientIter newIt = timeToCoefficient_.insert(it,std::pair<Time, KeyCoefficient>(time, KeyCoefficient(it->second.key, coefficient)));
  // Update keyToCoefficient_
  keyToCoefficient_[it->second.key] = newIt;
  // Remove the old coefficient
  timeToCoefficient_.erase(it);
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::removeCoefficientWithKey(Key key) {
  CHECK(hasCoefficientWithKey(key)) << "No coefficient with that key.";
  typename TimeToKeyCoefficientMap::iterator it1;
  typename boost::unordered_map<Key, CoefficientIter>::iterator it2;
  it2 = keyToCoefficient_.find(key);
  it1 = timeToCoefficient_.find(it1->second->first);
  timeToCoefficient_.erase(it1);
  keyToCoefficient_.erase(it2);
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::removeCoefficientAtTime(Time time) {
  CHECK(this->hasCoefficientAtTime(time)) << "No coefficient at that time.";
  typename TimeToKeyCoefficientMap::iterator it1;
  typename boost::unordered_map<Key, CoefficientIter>::iterator it2;
  it1 = timeToCoefficient_.find(time);
  it2 = keyToCoefficient_.find(it1->second.key);
  timeToCoefficient_.erase(it1);
  keyToCoefficient_.erase(it2);
}

/// \brief return true if there is a coefficient at this time
template <class Coefficient>
bool LocalSupport2CoefficientManager<Coefficient>::hasCoefficientAtTime(Time time) const {
  CoefficientIter it = timeToCoefficient_.find(time);
  return it != timeToCoefficient_.end();
}

/// \brief return true if there is a coefficient with this key
template <class Coefficient>
bool LocalSupport2CoefficientManager<Coefficient>::hasCoefficientWithKey(Key key) const {
  CoefficientIter it = keyToCoefficient_.find(key);
  return it != keyToCoefficient_.end();
}

/// \brief set the coefficient associated with this key
///
/// This function fails if there is no coefficient associated
/// with this key.
template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::updateCoefficientByKey(Key key, const Coefficient& coefficient) {
  typename  boost::unordered_map<Key, CoefficientIter>::iterator it = keyToCoefficient_.find(key);
  CHECK(it != keyToCoefficient_.end()) << "Key " << key << " is not in the container.";
  *const_cast<CoefficientType*>(&(it)->second->second.coefficient) = coefficient;
}

/// \brief get the coefficient associated with this key
template <class Coefficient>
Coefficient LocalSupport2CoefficientManager<Coefficient>::getCoefficientByKey(Key key) const {
  typename  boost::unordered_map<Key, CoefficientIter>::const_iterator it = keyToCoefficient_.find(key);
  CHECK(it != keyToCoefficient_.end() ) << "Key " << key << " is not in the container.";
  return it->second->second.coefficient;
}
template <class Coefficient>
Time LocalSupport2CoefficientManager<Coefficient>::getCoefficientTimeByKey(Key key) const {
  typename  boost::unordered_map<Key, CoefficientIter>::const_iterator it = keyToCoefficient_.find(key);
  CHECK(it != keyToCoefficient_.end()) << "Key " << key << " is not in the container.";
  return it->second->first;
}


/// \brief Get the coefficients that are active at a certain time.
template <class Coefficient>
bool LocalSupport2CoefficientManager<Coefficient>::getCoefficientsAt(Time time,
                                                                     CoefficientIter* outCoefficient0,
                                                                     CoefficientIter* outCoefficient1) const {
  CHECK_NOTNULL(outCoefficient0);
  CHECK_NOTNULL(outCoefficient1);
  if( timeToCoefficient_.empty() ) {
    LOG(INFO) << "No coefficients";
    return false;
  }

  CoefficientIter it;

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
  *outCoefficient0 = it;
  *outCoefficient1 = (++it);

  return true;
}

/// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::getCoefficientsInRange(
    Time startTime, Time endTime, CoefficientMap* outCoefficients) const {

  if (startTime <= endTime && startTime <= this->getMaxTime()
      && endTime >= this->getMinTime()) {
    // be forgiving if start time is lower than definition of curve
    if (startTime < this->getMinTime()) {
      startTime = this->getMinTime();
    }
    // be forgiving if end time is greater than definition of curve
    if (endTime > this->getMaxTime()) {
      endTime = this->getMaxTime();
    }
    CoefficientIter it;
    // set iterator to coefficient left or equal of start time
    it = timeToCoefficient_.upper_bound(startTime);
    it--;
    // iterate through coefficients
    for (; it != timeToCoefficient_.end() && it->first < endTime; ++it) {
      (*outCoefficients)[it->second->key] = it->second->coefficient;
    }
    if (it != timeToCoefficient_.end()) {
      (*outCoefficients)[it->second->key] = it->second->coefficient;
    }
  }
}

/// \brief Get all of the curve's coefficients.
template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::getCoefficients(CoefficientMap* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  CoefficientIter it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    (*outCoefficients)[it->second->key] = it->second->coefficient;
  }
}

/// \brief Set coefficients.
///
/// If any of these coefficients doen't exist, there is an error
template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::updateCoefficients(
    const CoefficientMap& coefficients) {
  typename CoefficientMap::const_iterator it;
  it = coefficients.cbegin();
  for (; it != coefficients.end(); ++it) {
    this->updateCoefficientByKey(it->first, it->second);
  }
}

/// \brief return the number of coefficients
template <class Coefficient>
Key LocalSupport2CoefficientManager<Coefficient>::size() const {
  return timeToCoefficient_.size();
}

template <class Coefficient>
bool LocalSupport2CoefficientManager<Coefficient>::empty() const {
  return timeToCoefficient_.empty();
}

/// \brief clear the coefficients
template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::clear() {
  keyToCoefficient_.clear();
  timeToCoefficient_.clear();
}

template <class Coefficient>
Time LocalSupport2CoefficientManager<Coefficient>::getMinTime() const {
  if (timeToCoefficient_.empty()) {
    return 0;
  }
  return timeToCoefficient_.begin()->first;
}

template <class Coefficient>
Time LocalSupport2CoefficientManager<Coefficient>::getMaxTime() const {
  if (timeToCoefficient_.empty()) {
    return 0;
  }
  return timeToCoefficient_.rbegin()->first;
}

template <class Coefficient>
void LocalSupport2CoefficientManager<Coefficient>::checkInternalConsistency(bool doExit) const {
  CHECK_EQ(keyToCoefficient_.size(), timeToCoefficient_.size());
  CoefficientIter it;
  typename boost::unordered_map<Key, CoefficientIter>::const_iterator itc;
  for(it = timeToCoefficient_.begin() ; it != timeToCoefficient_.end(); ++it) {
    Key key = it->second.key;
    itc = keyToCoefficient_.find(key);
    CHECK_EQ(itc->first, itc->second->second.key);
    CHECK( itc != keyToCoefficient_.end() ) << "Key " << key << " is not in the map";
    // This is probably the important one.
    // Check that the itc->second iterator
    // points to the same object as it->second.
    // It is supposedly guaranteed by the
    // unordered map interface that these
    // pointers never get reallocated.
    //todo the implementation differs from the comment above. Here only the
    //equality between coefficients is checked
    CHECK_EQ(itc->second->second.coefficient, it->second.coefficient);
  }
  if (doExit) {
    exit(0);
  }
}

template <class Coefficient>
bool LocalSupport2CoefficientManager<Coefficient>::hasCoefficientAtTime(Time time, CoefficientIter *it, double tol) {
  for ((*it) = timeToCoefficient_.begin();
      (*it) != timeToCoefficient_.end(); ++(*it)) {
    if ((*it)->first >= time-tol) {
      if ((*it)->first <= time+tol) {
        return true;
      } else {
        return false;
      }
    }
  }
  return false;
}

} // namespace
