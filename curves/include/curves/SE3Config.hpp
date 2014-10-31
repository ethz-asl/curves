#ifndef SE3CONFIG_H_
#define SE3CONFIG_H_

#include <Eigen/Core>
#include "kindr/minimal/quat-transformation.h"

namespace curves {

typedef Eigen::Matrix<double, 6, 1> Vector6d;

struct SE3Config {
  typedef kindr::minimal::QuatTransformationTemplate<double> ValueType;
  typedef Vector6d DerivativeType;
};

}  // namespace curves

#endif // SE3CONFIG_H_
