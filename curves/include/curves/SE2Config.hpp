/*
 * @file SE2Config.hpp
 * @date Nov 24, 2015
 * @author Renaud Dubé
 */

#ifndef SE2CONFIG_H_
#define SE2CONFIG_H_

#include <Eigen/Core>
#include "gtsam/geometry/Pose2.h"

namespace curves {

typedef Eigen::Matrix<double, 3, 1> Vector3d;

struct SE2Config {
  typedef gtsam::Pose2 ValueType;
  typedef Vector3d DerivativeType;
};

}  // namespace curves

#endif // SE2CONFIG_H_
