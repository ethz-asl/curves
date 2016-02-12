/*
 * @file SE3CurveFactory.hpp
 * @date Feb 4th, 2016
 * @author Renaud Dubé, Abel Gawel
 */

#ifndef SE3_CURVE_FACTORY_HPP_
#define SE3_CURVE_FACTORY_HPP_

#include "SE3Curve.hpp"
#include "SlerpSE3Curve.hpp"
#include "DiscreteSE3Curve.hpp"
#include "SemiDiscreteSE3Curve.hpp"
#include "SE3CompositionCurve.hpp"

#include <string>

namespace curves {

class SE3CurveFactory {

 public:

  SE3CurveFactory() {};
  ~SE3CurveFactory() {};

  static std::shared_ptr<SE3Curve> create_curve(const std::string& curveType);

}; // class SE3CurveFactory
} // namespace ctsm

#endif // SE3_CURVE_FACTORY_HPP_
