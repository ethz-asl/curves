#ifndef CURVES_CURVE_HPP
#define CURVES_CURVE_HPP

#include "types.hpp"
#include "Coefficient.hpp"

namespace curves {

class CurveBase {
 public:
  
  CurveBase();
  virtual ~CurveBase();

  ///\defgroup Info 
  ///\name Methods to get information about the curve.
  ///@{

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const = 0;

  /// \brief The dimension of the underlying manifold
  virtual size_t dim() const = 0;

  ///@}

  /// 
  /// \defgroup GetSet 
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


  /// \defgroup Time 
  /// \name Methods to deal with time.
  ///@{

  /// The first valid time of the curve.
  virtual Time getMinTime() const = 0;
  
  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const = 0;

  ///@}

  /// \defgroup GrowShrink 
  /// \name Methods to grow or shrink the curve.
  ///@{

  /// \brief extend or retract the curve to the specified range, \f$[t_s,t_e)\f$.
  virtual void setTimeRange(Time minTime, Time maxTime);
  
  ///@}


};

} // namespace curves



#endif /* CURVES_CURVE_HPP */
