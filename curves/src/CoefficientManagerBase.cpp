/*
 * @file CoefficientManagerBase.cpp
 * @date Dec 3, 2014
 * @author Paul Furgale, Abel Gawel, Renaud Dube
 */

#include <curves/CoefficientManagerBase.hpp>
#include <curves/KeyGenerator.hpp>
#include <boost/make_shared.hpp>
#include <glog/logging.h>

namespace curves {

CoefficientManagerBase::CoefficientManagerBase() {
}

CoefficientManagerBase::CoefficientManagerBase(const CoefficientManagerBase& rhs) {
  // Iterate over sorted rhs map
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator rhsIter = rhs.timeToCoefficient_.begin();
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::iterator hintIter = this->timeToCoefficient_.begin();
  for( ; rhsIter != rhs.timeToCoefficient_.end(); ++rhsIter) {
    // Deep copy the shared_ptr and clone() for possible inheritences of KeyCoefficientTime
    boost::shared_ptr<KeyCoefficientTime> temp = boost::shared_ptr<KeyCoefficientTime>(rhsIter->second->clone());
    CHECK_EQ(rhsIter->first, temp->time);

    // Insert into maps
    hintIter = this->timeToCoefficient_.insert(hintIter, std::make_pair(rhsIter->first, temp));
    typedef boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::iterator iterator;
    std::pair<iterator, bool> success = this->coefficients_.insert(std::make_pair(temp->key, temp));
    CHECK(success.second);
  }
}

CoefficientManagerBase& CoefficientManagerBase::operator= (CoefficientManagerBase rhs) {
    std::swap(this->coefficients_,      rhs.coefficients_);
    std::swap(this->timeToCoefficient_, rhs.timeToCoefficient_);
    return *this;
}

CoefficientManagerBase::~CoefficientManagerBase() {
}

/// Compare this Coeficient with another for equality.
bool CoefficientManagerBase::equals(const CoefficientManagerBase& other,
                                    double tol) const {
  bool equal = true;
  equal &= coefficients_.size() == other.coefficients_.size();
  equal &= timeToCoefficient_.size() == other.timeToCoefficient_.size();
  if (equal) {
    boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it1, it2;
    it1 = coefficients_.begin();
    it2 = other.coefficients_.begin();
    for( ; it1 != coefficients_.end(); ++it1, ++it2) {
      equal &= it1->first == it2->first;
      if (it1->second && it2->second) {
        equal &= it1->second->equals(*(it2->second));
      } else {
        equal &= it1->second == it2->second;
      }
    }
    std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it3, it4;
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

void CoefficientManagerBase::getKeys(std::vector<Key>* outKeys) const {
  CHECK_NOTNULL(outKeys);
  outKeys->clear();
  appendKeys(outKeys);
}

void CoefficientManagerBase::appendKeys(std::vector<Key>* outKeys) const {
  CHECK_NOTNULL(outKeys);
  outKeys->reserve(outKeys->size() + coefficients_.size());
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    outKeys->push_back(it->second->key);
  }
}

void CoefficientManagerBase::getTimes(std::vector<Time>* outTimes) const {
  CHECK_NOTNULL(outTimes);
  outTimes->clear();
  outTimes->reserve(timeToCoefficient_.size());
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    outTimes->push_back(it->first);
  }
}

void CoefficientManagerBase::print(const std::string& str) const {
  // \todo (Abel or Renaud)
}

Key CoefficientManagerBase::insertCoefficient(Time time, const Coefficient& coefficient) {
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::iterator it;
  Key key;
  if (this->hasCoefficientAtTime(time, &it)) {
    this->setCoefficientByKey(it->second->key, coefficient);
    key = it->second->key;
  } else {
    key = KeyGenerator::getNextKey();
    this->insertNewCoefficient(key, time, coefficient);
  }
  return key;
}

void CoefficientManagerBase::insertNewCoefficient(Key key, Time time, const Coefficient& coefficient) {
  typedef boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::iterator iterator;
  std::pair<iterator, bool> success = coefficients_.insert(
      std::make_pair(key, makeContainer(key, coefficient, time)));
  timeToCoefficient_[time] = success.first->second;
}

boost::shared_ptr<KeyCoefficientTime> CoefficientManagerBase::makeContainer(Key key, Coefficient coefficient, Time time) {
  return boost::make_shared<KeyCoefficientTime>(key, coefficient, time);
}

/// \brief insert coefficients. Optionally returns the keys for these coefficients
void CoefficientManagerBase::insertCoefficients(const std::vector<Time>& times,
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

/// \brief return true if there is a coefficient at this time
bool CoefficientManagerBase::hasCoefficientAtTime(Time time) const {
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it = timeToCoefficient_.find(time);
  return it != timeToCoefficient_.end();
}

/// \brief return true if there is a coefficient with this key
bool CoefficientManagerBase::hasCoefficientWithKey(Key key) const {
  boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it = coefficients_.find(key);
  return it != coefficients_.end();
}

/// \brief set the coefficient associated with this key
///
/// This function fails if there is no coefficient associated
/// with this key.
void CoefficientManagerBase::setCoefficientByKey(Key key, const Coefficient& coefficient) {
  boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::iterator it = coefficients_.find(key);
  CHECK( it != coefficients_.end() ) << "Key " << key << " is not in the container.";
  it->second->coefficient = coefficient;
}

/// \brief set the coefficient associated with this key
///
/// This function fails if there is no coefficient associated
/// with this key.
void CoefficientManagerBase::setCoefficientVectorByKey(Key key, const Eigen::VectorXd& vector) {
  boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::iterator it = coefficients_.find(key);
  CHECK( it != coefficients_.end() ) << "Key " << key << " is not in the container.";
  it->second->coefficient.setVector(vector);
}

/// \brief get the coefficient associated with this key
Coefficient CoefficientManagerBase::getCoefficientByKey(Key key) const {
  boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it = coefficients_.find(key);
  CHECK( it != coefficients_.end() ) << "Key " << key << " is not in the container.";
  return it->second->coefficient;
}

/// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
void CoefficientManagerBase::getCoefficientsInRange(
    Time startTime, Time endTime, Coefficient::Map* outCoefficients) const {

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
    std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it;
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
void CoefficientManagerBase::getCoefficients(Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    (*outCoefficients)[it->second->key] = it->second->coefficient;
  }
}

/// \brief Set coefficients.
///
/// If any of these coefficients doen't exist, there is an error
void CoefficientManagerBase::setCoefficients(
    const Coefficient::Map& coefficients) {
  boost::unordered_map<Key, Coefficient>::const_iterator it;
  it = coefficients.begin();
  for (; it != coefficients.end(); ++it) {
    this->setCoefficientByKey(it->first, it->second);
  }
}

/// \brief return the number of coefficients
Key CoefficientManagerBase::size() const {
  return coefficients_.size();
}

/// \brief clear the coefficients
void CoefficientManagerBase::clear() {
  coefficients_.clear();
  timeToCoefficient_.clear();
}

Time CoefficientManagerBase::getMinTime() const {
  if (timeToCoefficient_.empty()) {
    return 0;
  }
  return timeToCoefficient_.begin()->first;
}

Time CoefficientManagerBase::getMaxTime() const {
  if (timeToCoefficient_.empty()) {
    return 0;
  }
  return timeToCoefficient_.rbegin()->first;
}

void CoefficientManagerBase::checkInternalConsistency(bool doExit) const {
  CHECK_EQ(coefficients_.size(), timeToCoefficient_.size());
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::const_iterator it;
  it = timeToCoefficient_.begin();
  for( ; it != timeToCoefficient_.end(); ++it) {
    CHECK(it->second);
    CHECK_EQ(it->first, it->second->time);
    Key key = it->second->key;
    boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> >::const_iterator itc = coefficients_.find(key);
    CHECK( itc != coefficients_.end() ) << "Key " << key << " is not in the map";
    // This is probably the important one.
    // Check that the it->second pointer
    // points to the same object as itc->second.
    // It is supposedly guaranteed by the 
    // unordered map interface that these
    // pointers never get reallocated.
    CHECK_EQ(itc->second, it->second);
  }
  if (doExit) {
    exit(0);
  }
}

void CoefficientManagerBase::removeCoefficientWithKey(Key key) {
  CHECK(false) << "Not implemented";
  // todo Abel and Renaud
}

void CoefficientManagerBase::removeCoefficientAtTime(Time time) {
  CHECK(false) << "Not implemented";
  // todo Abel and Renaud
}

bool CoefficientManagerBase::hasCoefficientAtTime(Time time, std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::iterator *it) {
  *it = timeToCoefficient_.find(time);
  return *it != timeToCoefficient_.end();
}

} // namespace curves
