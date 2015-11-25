/*
 * @file Pose2_Expressions.hpp
 * @date Nov 24, 2015
 * @author Renaud Dubé, Abel Gawel
 */

#ifndef POSE2_EXPRESSIONS_HPP
#define POSE2_EXPRESSIONS_HPP

#include "gtsam/nonlinear/Expression.h"
#include "gtsam/geometry/Pose2.h"

typedef gtsam::Pose2 SE2;
typedef gtsam::Expression<SE2> ETransformation;
typedef gtsam::OptionalJacobian<3, 3> ChartJacobian;
typedef gtsam::Vector3 Vector3;
typedef gtsam::Expression<gtsam::Vector3> EVector3;

SE2 inverseImplementation(
    const SE2& T, ChartJacobian HT) {
  return T.inverse(HT);
}

ETransformation inverse(
    const ETransformation& T) {
  return ETransformation(&inverseImplementation, T);
}

SE2 composeImplementation(const SE2& T1, const SE2& T2,
                          ChartJacobian HT1, ChartJacobian HT2) {
  return T1.compose(T2, HT1, HT2);
}

ETransformation compose(const ETransformation& T1, const ETransformation& T2) {
  return ETransformation(&composeImplementation, T1, T2);
}

Vector3 transformationLogImplementation(const SE2& T, ChartJacobian HT) {
  return gtsam::Pose2::Logmap(T, HT);
}

EVector3 transformationLog(const ETransformation& T) {
  return EVector3(&transformationLogImplementation, T);
}

SE2 transformationExpImplementation(const Vector3& params,
                                                ChartJacobian Hp) {
  return gtsam::Pose2::Expmap(params, Hp);
}

ETransformation transformationExp(const EVector3& params) {
  return ETransformation(&transformationExpImplementation, params);
}

template <int N>
Eigen::Matrix<double, N, 1> vectorScalingImplementation(const Eigen::Matrix<double, N, 1> & v, double alpha,
                                                        gtsam::OptionalJacobian<N, N> H1,
                                                        gtsam::OptionalJacobian<N, 1> H2) {
  if (H1) {
    *H1 = gtsam::OptionalJacobian<N,N>::Jacobian::Identity()*alpha;
  }
  if (H2) {
    *H2 = v;
  }
  return v*alpha;
}

template <int N>
gtsam::Expression<Eigen::Matrix<double, N, 1> > vectorScaling(const gtsam::Expression<Eigen::Matrix<double, N, 1> >& v, double alpha) {
  return gtsam::Expression<Eigen::Matrix<double, N, 1> >(boost::bind(&vectorScalingImplementation<N>, _1, alpha, _2, boost::none), v);
}

ETransformation slerp(
    const ETransformation& T0,
    const ETransformation& T1,
    double alpha) {
  return compose(T0, transformationExp(vectorScaling(transformationLog(compose(inverse(T0), T1)), alpha)));
}

#endif /* POSE2_EXPRESSIONS_HPP */
