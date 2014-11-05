#ifndef SE3CONFIG_H_
#define SE3CONFIG_H_

#include <Eigen/Core>
#include "kindr/minimal/quat-transformation.h"

#ifdef CURVES_USE_GTSAM
#include <gtsam/base/Manifold.h>
#endif

namespace curves {

typedef Eigen::Matrix<double, 6, 1> Vector6d;

struct SE3Config {
  typedef kindr::minimal::QuatTransformationTemplate<double> ValueType;
  typedef Vector6d DerivativeType;
};

}  // namespace curves

#ifdef CURVES_USE_GTSAM

namespace gtsam {
namespace traits {

using curves::SE3Config;

template<>
struct is_group<SE3Config::ValueType> : public boost::true_type {
};

template<>
struct is_manifold<SE3Config::ValueType> : public boost::true_type {
};

template<>
struct dimension<SE3Config::ValueType> : public boost::integral_constant<int, 6> {
};

}
}

#endif


#endif // SE3CONFIG_H_
