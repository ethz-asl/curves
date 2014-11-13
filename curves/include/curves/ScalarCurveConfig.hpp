/*
 * @file ScalarCurveConfig.hpp
 * @date Aug 18, 2014
 * @author Paul Furgale
 */

#ifndef CURVES_SCALAR_CURVE_CONFIG_HPP
#define CURVES_SCALAR_CURVE_CONFIG_HPP

#include <Eigen/Core>

namespace curves {

struct ScalarCurveConfig {
  typedef double ValueType;
  typedef double DerivativeType;
};

} // namespace curves


#endif /* CURVES_SCALAR_CURVE_CONFIG_HPP */
