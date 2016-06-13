/*
 * @file Pose2_Expressions.hpp
 * @date Nov 24, 2015
 * @author Renaud Dubé, Abel Gawel
 *
 * This helper implements expressions for Pose2 objects in order to
 * perform spherical linear interpolation (slerp) in SE2.
 */

#ifndef POSE2_EXPRESSIONS_HPP
#define POSE2_EXPRESSIONS_HPP

#include "gtsam/nonlinear/Expression.h"
#include "gtsam/geometry/Pose2.h"
#include "gtsam/geometry/Point2.h"

#define N_TEST_ITERATIONS 10000

typedef gtsam::Pose2 SE2;
typedef gtsam::Expression<SE2> ETransformation;
typedef gtsam::OptionalJacobian<3, 3> ChartJacobian;
typedef gtsam::Vector3 Vector3;
typedef gtsam::Expression<gtsam::Vector3> EVector3;

inline SE2 makeRandomSE2() {
  gtsam::Vector2 pval;
  pval.setRandom();

  return SE2(((double)rand() / RAND_MAX)*6.2832, gtsam::Point2(pval));
}

inline SE2 inverseImplementation(
    const SE2& T, ChartJacobian HT) {
  return T.inverse(HT);
}

inline ETransformation inverse(
    const ETransformation& T) {
  return ETransformation(&inverseImplementation, T);
}

inline SE2 composeImplementation(const SE2& T1, const SE2& T2,
                                 ChartJacobian HT1, ChartJacobian HT2) {
  return T1.compose(T2, HT1, HT2);
}

inline ETransformation compose(const ETransformation& T1, const ETransformation& T2) {
  return ETransformation(&composeImplementation, T1, T2);
}

inline double function1(double alpha) {
  double K = 0.0000330688 / 0.0833333;
  double threshold = pow((std::numeric_limits<double>::epsilon()/K),0.25)*10;
  double rval;
  if (fabs(alpha) > threshold) {
    rval = 1/alpha-0.5*sin(alpha)/(1-cos(alpha));
  } else {
    rval = 0.0833333 * alpha + 0.00138889 * pow(alpha,3);
  }
  return rval;
}

inline double function2(double alpha) {
  double K = 0.00138889;
  double threshold = pow((std::numeric_limits<double>::epsilon()/K),0.25)*10;
  double rval;
  if (fabs(alpha) > threshold) {
    rval = alpha*0.5*sin(alpha)/(1-cos(alpha));
  } else {
    rval = 1 - 0.0833333 * pow(alpha,2);
  }
  return rval;
}

inline Vector3 transformationLogImplementation(const SE2& T, ChartJacobian HT) {
  Vector3 v = gtsam::Pose2::Logmap(T, HT);
  if (HT) {
    double alpha = v[2];
    gtsam::Matrix3 J;

    double alphaInv = 1/alpha;
    double halfCotHalfAlpha = 0.5*sin(alpha)/(1-cos(alpha));
    double v1 = v[0], v2 = v[1];
    J << function2(alpha), -0.5*alpha, v1*function1(alpha) + 0.5*v2,
        0.5*alpha, function2(alpha),  v2*function1(alpha) - 0.5*v1,
        0, 0, 1;

    *HT = J;
  }
  return v;
}

inline EVector3 transformationLog(const ETransformation& T) {
  return EVector3(&transformationLogImplementation, T);
}

inline SE2 transformationExpImplementation(const Vector3& params,
                                           ChartJacobian Hp) {
  return gtsam::Pose2::Expmap(params, Hp);
}

inline ETransformation transformationExp(const EVector3& params) {
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

inline ETransformation slerp(
    const ETransformation& T0,
    const ETransformation& T1,
    double alpha) {
  return compose(T0, transformationExp(vectorScaling(transformationLog(compose(inverse(T0), T1)), alpha)));
}

inline gtsam::Point2 transformFromCurveImplementation(gtsam::Pose2 pose, gtsam::Point2 point,
                                                      gtsam::OptionalJacobian<2, 3> H1,
                                                      gtsam::OptionalJacobian<2, 2> H2) {
  return pose.transform_from(point, H1, H2);
}

inline gtsam::Expression<gtsam::Point2> transformFromCurve(const gtsam::Expression<gtsam::Pose2>& Epose,
                                                           const gtsam::Expression<gtsam::Point2>& Epoint) {
  return gtsam::Expression<gtsam::Point2>(&transformFromCurveImplementation,
                                          Epose, Epoint);
}

inline double distanceBetweenPointsImplementation(gtsam::Point2 p1, gtsam::Point2 p2,
                                                  gtsam::OptionalJacobian<1, 2> H1,
                                                  gtsam::OptionalJacobian<1, 2> H2) {
  return p1.distance(p2, H1, H2);
}

inline gtsam::Expression<double> distanceBetweenPoints(const gtsam::Expression<gtsam::Point2>& Ep1,
                                                       const gtsam::Expression<gtsam::Point2>& Ep2) {
  return gtsam::Expression<double>(&distanceBetweenPointsImplementation,
                                   Ep1, Ep2);
}

inline gtsam::Point2 pointsSubtractionImplementation(gtsam::Point2 p1, gtsam::Point2 p2,
                                                     gtsam::OptionalJacobian<2, 2> H1,
                                                     gtsam::OptionalJacobian<2, 2> H2) {
  if (H1) {
    *H1 << 1, 0, 0, 1;
  }

  if (H2) {
    *H2 << -1, 0, 0, -1;
  }

  return gtsam::Point2(p1.x()-p2.x(),p1.y()-p2.y());
}

inline gtsam::Expression<gtsam::Point2> pointsSubtraction(const gtsam::Expression<gtsam::Point2>& Ep1,
                                                          const gtsam::Expression<gtsam::Point2>& Ep2) {
  return gtsam::Expression<gtsam::Point2>(&pointsSubtractionImplementation,
                                          Ep1, Ep2);
}

inline gtsam::Point2 pose2AsPoint2Implementation(gtsam::Pose2 pose, gtsam::OptionalJacobian<2, 3> H) {
  if (H) {
    *H << 1, 0, 0, 0, 1, 0;
  }
  return gtsam::Point2(pose.x(), pose.y());
}

inline gtsam::Expression<gtsam::Point2> pose2AsPoint2(const gtsam::Expression<gtsam::Pose2>& Epose) {
  return gtsam::Expression<gtsam::Point2>(&pose2AsPoint2Implementation, Epose);
}

inline double poseRangeImplementation(const SE2& T1, const SE2& T2,
                                  gtsam::OptionalJacobian<1, 3> H1,
                                  gtsam::OptionalJacobian<1, 3> H2) {
  return T1.range(T2, H1, H2);
}

inline gtsam::Expression<double> poseRange(const gtsam::Expression<gtsam::Pose2>& Epose1,
                                           const gtsam::Expression<gtsam::Pose2>& Epose2) {
  return gtsam::Expression<double>(&poseRangeImplementation, Epose1, Epose2);
}

inline double pointRangeImplementation(const SE2& pose, const gtsam::Point2& point,
                                  gtsam::OptionalJacobian<1, 3> H1,
                                  gtsam::OptionalJacobian<1, 2> H2) {
  return pose.range(point, H1, H2);
}

inline gtsam::Expression<double> pointRange(const gtsam::Expression<gtsam::Pose2>& Epose,
                                           const gtsam::Expression<gtsam::Point2>& Epoint) {
  return gtsam::Expression<double>(&pointRangeImplementation, Epose, Epoint);
}

#endif /* POSE2_EXPRESSIONS_HPP */
