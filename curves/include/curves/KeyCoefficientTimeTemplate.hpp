/*
 * @file KeyCoefficientTime.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale
 */

#ifndef CT_KEY_COEFFICIENT_HPP
#define CT_KEY_COEFFICIENT_HPP

#include <boost/cstdint.hpp>

namespace curves {

typedef boost::int64_t Time;   // RD Why was Types.hpp removed?
typedef size_t Key;

/// \struct KeyCoefficientTime
/// \brief Stores a key, coefficient, and time all together
///
/// A helper struct for hermite-type curves
template <class Coefficient>
struct KeyCoefficientTimeTemplate {
  typedef Coefficient CoefficientType;

  Key key;
  Coefficient coefficient;
  Time time;
  
  KeyCoefficientTimeTemplate(Key key, const Coefficient& coefficient, Time time) :
      key(key), coefficient(coefficient), time(time) {}
  KeyCoefficientTimeTemplate() {};
  bool equals(const KeyCoefficientTimeTemplate& other) const {
    return key == other.key && time == other.time &&
        coefficient.equals(other.coefficient);
  }
};

} // namespace curves


#endif /* CT_KEY_COEFFICIENT_HPP */
