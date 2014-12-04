/*
 * @file KeyCoefficientTime.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale
 */

#ifndef CT_KEY_COEFFICIENT_HPP
#define CT_KEY_COEFFICIENT_HPP

#include "types.hpp"
#include "Coefficient.hpp"

namespace curves {

/// \struct KeyCoefficientTime
/// \brief Stores a key, coefficient, and time all together
///
/// A helper struct for hermite-type curves
struct KeyCoefficientTime {
  Key key;
  Coefficient coefficient;
  Time time;
  
  KeyCoefficientTime(Key key, const Coefficient& coefficient, Time time) :
      key(key), coefficient(coefficient), time(time) {}
  KeyCoefficientTime() {};
  virtual bool equals(const KeyCoefficientTime& other) const {
    return key == other.key && time == other.time &&
        coefficient.equals(other.coefficient);
  }
};

} // namespace curves


#endif /* CT_KEY_COEFFICIENT_HPP */
