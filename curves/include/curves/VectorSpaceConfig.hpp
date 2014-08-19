#ifndef CURVES_VECTOR_SPACE_CONFIG_HPP
#define CURVES_VECTOR_SPACE_CONFIG_HPP

#include <Eigen/Core>

namespace curves {

struct VectorSpaceConfig {
  typedef Eigen::VectorXd ValueType;
  typedef Eigen::VectorXd DerivativeType;
};

} // namespace curves


#endif /* CURVES_VECTOR_SPACE_CONFIG_HPP */
