#ifndef CURVES_CURVE_HPP
#define CURVES_CURVE_HPP

#include "types.hpp"
#include "Coefficient.hpp"

namespace curves {

class Curve {
 public:
  
  Curve();
  virtual ~Curve();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const = 0;
  
  /// \name Methods to get and set coefficients.
  ///@{

  /// \brief Get the coefficients that are active at a certain time.
  virtual void getCoefficientsAt(Time time, 
                                 Coefficient::Map& outCoefficients) const = 0;

  /// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
  virtual void getCoefficientsInRange(Time startTime, 
                                      Time endTime, 
                                      Coefficient::Map& outCoefficients) const = 0;

  /// \brief Get all of the curve's coefficients.
  virtual void getCoefficients(Coefficient::Map& outCoefficients) const = 0;
  
  /// \brief Set a coefficient.
  virtual void setCoefficient(Key key, const Coefficient& value) = 0;

  /// \brief Set coefficients.
  virtual void setCoefficients(Coefficient::Map& coefficients) = 0;
  ///@}


  /// \name Methods to deal with time.
  ///@{

  /// The first valid time for the curve.
  virtual Time getBackTime() const = 0;
  
  /// The one past the last valid time for the curve.
  virtual Time getFrontTime() const = 0;

  ///@}

  /// \name Methods to grow or shrink the curve.
  ///@{

  /// Extend the curve into the future so that it can be evaluated
  /// at this time.
  virtual void extendFront(Time time) = 0;

  /// Extend the curve into the past so that it can be evaluated
  /// at this time.
  virtual void extendBack(Time time) = 0;

  /// Retract the curve in the front. The caller guarantees that
  /// they will never try to evaluate this time or higher again.
  virtual void retractFront(Time time) = 0;

  /// Retract the curve in the back. The caller guarantees that
  /// they will never try to evaluate this time or lower again.
  virtual void retractBack(Time time) = 0;

  ///@}

  
  
};

} // namespace curves



#endif /* CURVES_CURVE_HPP */
