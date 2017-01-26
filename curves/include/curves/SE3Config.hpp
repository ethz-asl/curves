/*
 * ScalarCurveConfig.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Paul Furgale, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include <Eigen/Core>
#include <kindr/Core>

namespace curves {

typedef Eigen::Matrix<double, 6, 1> Vector6d;

struct SE3Config {
  typedef kindr::HomTransformQuatD ValueType;
  typedef kindr::TwistGlobalD DerivativeType;
};

} // namespace
