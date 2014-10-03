#ifndef SE3CONFIG_H_
#define SE3CONFIG_H_

namespace curves {

typedef Eigen::Matrix<double, 6, 1> Vector6d;

struct SE3Config {
  typedef Eigen::Matrix4d ValueType;
  typedef Vector6d DerivativeType;
};

}  // namespace curves

#endif // SE3CONFIG_H_
