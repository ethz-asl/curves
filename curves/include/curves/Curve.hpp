/*
 * @file Curve.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale, Renaud Dube
 */

#ifndef CURVES_CURVE_HPP
#define CURVES_CURVE_HPP

#include "gtsam/nonlinear/Expression.h"
#include <boost/cstdint.hpp>

namespace curves {

typedef boost::int64_t Time;
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
  virtual ValueType evaluate(Time time) const = 0;
  
  /// Evaluate the curve derivatives.
  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) const = 0;

//  /// \brief Get an evaluator at this time.
//  virtual EvaluatorTypePtr getEvaluator(const Time& time) const = 0;

  /// \brief Get a gtsam::Expression which evaluates the curve at this time.
  virtual gtsam::Expression<ValueType> getValueExpression(const Time& time) const = 0;

  /// \brief Get a gtsam::Expression which evaluates the derivative of the curve at this time.
  virtual gtsam::Expression<DerivativeType> getDerivativeExpression(const Time& time, unsigned derivativeOrder) const = 0;

  ///@}

  /// \name Methods to fit the curve based on data.
  ///@{

  /// Extend the curve so that it can be evaluated at these times.
  /// Try to make the curve fit to the values.
  /// Underneath the curve should have some default policy for fitting.
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys) = 0;

  /// \brief Fit a new curve to these data points.
  ///
  /// The existing curve will be cleared.
  /// Underneath the curve should have some default policy for fitting.
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL) = 0;

  ///@}

   /// Initialize a GTSAM values structure with the desired keys
   virtual void initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) const = 0;

   /// Initialize a GTSAM values structure for all keys
   virtual void initializeGTSAMValues(gtsam::Values* values) const = 0;

   // updates the relevant curve coefficients from the GTSAM values structure
   virtual void updateFromGTSAMValues(const gtsam::Values& values) = 0;

   virtual void clear() = 0;

   /// \brief Perform a rigid transformation on the left side of the curve
   virtual void transformCurve(const ValueType T) = 0;

};

} // namespace curves


#endif /* CURVES_CURVE_HPP */
