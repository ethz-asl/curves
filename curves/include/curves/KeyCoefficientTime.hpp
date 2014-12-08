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
  KeyCoefficientTime() {}
  virtual KeyCoefficientTime* clone() const {
    return new KeyCoefficientTime(*this);
  }
  virtual bool equals(const KeyCoefficientTime& other) const {
    return key == other.key && time == other.time &&
        coefficient.equals(other.coefficient);
  }
};

/// \struct KeyCoefficientTime_CRTP
/// \brief This inheritence-helper class implements clone() for all derived classes
///
/// CRTP class for cloning (see http://en.wikipedia.org/wiki/Curiously_recurring_template_pattern)
template <typename Derived>
struct KeyCoefficientTime_CRTP : public KeyCoefficientTime {

  KeyCoefficientTime_CRTP(Key key, const Coefficient& coefficient, Time time) :
    KeyCoefficientTime(key, coefficient, time) {}
  KeyCoefficientTime_CRTP() {}

  virtual KeyCoefficientTime* clone() const {
    return new Derived(static_cast<Derived const&>(*this));
  }
};

} // namespace curves

#endif /* CT_KEY_COEFFICIENT_HPP */
