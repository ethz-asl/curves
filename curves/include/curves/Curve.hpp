/*
 * Curve.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Paul Furgale, Abel Gawel, Renaud Dube, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include <string>
#include <vector>

namespace curves {

typedef double Time;
typedef size_t Key;

template<typename CurveConfig>
class Curve
{
 public:

  /// The value type of the curve.
  typedef typename CurveConfig::ValueType ValueType;

  /// The curve's derivative type.
  typedef typename CurveConfig::DerivativeType DerivativeType;

  Curve() { }
  virtual ~Curve() { }

  ///\defgroup Info
  ///\name Methods to get information about the curve.
  ///@{

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const = 0;

  /// \brief The dimension of the underlying manifold
  //size_t dim() const;  // get this form the curve's value type

  /// The first valid time of the curve.
  virtual Time getMinTime() const = 0;

  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const = 0;
  ///@}

  /// \name Methods to evaluate the curve
  ///@{

  /// Evaluate the ambient space of the curve.
  virtual bool evaluate(ValueType& value, Time time) const = 0;

//  /// Evaluate the curve derivatives.
  virtual bool evaluateDerivative(DerivativeType& derivative, Time time, unsigned derivativeOrder) const = 0;

  ///@}

  /// \name Methods to fit the curve based on data.
  ///@{

  /// Extend the curve so that it can be evaluated at these times.
  /// Try to make the curve fit to the values.
  /// Underneath the curve should have some default policy for fitting.
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys = NULL) = 0;

  /// \brief Fit a new curve to these data points.
  ///
  /// The existing curve will be cleared.
  /// Underneath the curve should have some default policy for fitting.
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL) = 0;

  ///@}

   /// \brief Clear all the curve coefficients
   virtual void clear() = 0;

   /// \brief Perform a rigid transformation on the left side of the curve
   virtual void transformCurve(const ValueType T) = 0;
};

}  // namespace
