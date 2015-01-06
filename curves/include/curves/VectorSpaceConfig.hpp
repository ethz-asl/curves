/*
 * @file VectorSpaceConfig.hpp
 * @date Aug 18, 2014
 * @author Paul Furgale, Renaud Dube
 */

#ifndef CURVES_VECTOR_SPACE_CONFIG_HPP
#define CURVES_VECTOR_SPACE_CONFIG_HPP

#include <Eigen/Core>

namespace curves {

template <int N>
struct VectorSpaceConfig {
  typedef Eigen::Matrix<double,N,1> ValueType;
  typedef Eigen::Matrix<double,N,1> DerivativeType;
};

} // namespace curves


#endif /* CURVES_VECTOR_SPACE_CONFIG_HPP */
