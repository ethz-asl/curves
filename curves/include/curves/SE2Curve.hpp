/*
 * @file SE2Curve.hpp
 * @date Nov 24, 2015
 * @author Renaud Dubé
 */

#ifndef SE2_CURVE_H_
#define SE2_CURVE_H_

#include "SE2Config.hpp"
#include "Curve.hpp"

namespace curves {

// Curves over SE2 inherit the interface from Curve and CurveBase and define specific
// methods to support physical interpretations of temporal derivatives.
//
// For the purposes of these uses, the curve is defined between Frame a and Frame b
// such that evaluate() returns \f$ \mathbf{T}_{a,b} \f$, the transformation that takes
// points from Frame b to Frame a.
//
class SE2Curve : public Curve<SE2Config> {
 public:
  SE2Curve();
  virtual ~SE2Curve();

  typedef Curve<SE2Config> Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;

};

}  // namespace curves

#endif // SE2_CURVE_H_
