/*
 * @file CubicHermiteSE3Curve.hpp
 * @date Feb 10, 2015
 * @author Abel Gawel, Renaud Dube
 */

#ifndef CURVES_CUBIC_HERMITE_SE3_CURVE_HPP
#define CURVES_CUBIC_HERMITE_SE3_CURVE_HPP

#include "SE3Curve.hpp"
#include "LocalSupport2CoefficientManager.hpp"
#include "kindr/minimal/cubic-hermite-translation-gtsam.h"
#include "kindr/minimal/cubic-hermite-quaternion-gtsam.h"
#include "kindr/minimal/rotation-quaternion-gtsam.h"
#include "kindr/minimal/common-gtsam.h"

// wrapper class for Hermite-style coefficients (made of QuatTransformation and Vector6)
namespace kindr {
namespace minimal {

typedef gtsam::Expression<kindr::minimal::RotationQuaternion> EQuaternion;
typedef gtsam::Expression<Eigen::Vector3d> EVector3;
typedef gtsam::Expression<kindr::minimal::QuatTransformation> ETransformation;

template <typename Scalar>
struct HermiteTransformation {
  typedef QuatTransformationTemplate<Scalar> QuatTransformation;
  typedef Eigen::Matrix<Scalar, 6, 1> Vector6;

 public:
  HermiteTransformation();
  HermiteTransformation(const QuatTransformation& transform, const Vector6& derivatives);
  ~HermiteTransformation();

  QuatTransformation getTransformation() const {
    return transformation_;
  }

  Vector6 getTransformationDerivative() const {
    return transformationDerivative_;
  }

  void setTransformation(const QuatTransformation& transformation) {
    transformation_ = transformation;
  }

  void setTransformationDerivative(const Vector6& transformationDerivative) {
    transformationDerivative_ = transformationDerivative;
  }

 private:
  QuatTransformation transformation_;
  Vector6 transformationDerivative_;
};
template <typename Scalar>
HermiteTransformation<Scalar>::HermiteTransformation() {};

template <typename Scalar>
HermiteTransformation<Scalar>::HermiteTransformation(const QuatTransformation& transform,
                                                     const Vector6& derivatives) :
                                                     transformation_(transform),
                                                     transformationDerivative_(derivatives) {};

template <typename Scalar>
HermiteTransformation<Scalar>::~HermiteTransformation() {};

}
}

// traits for Hermite-style coefficients wrapper
namespace gtsam {
template<> struct traits<kindr::minimal::HermiteTransformation<double>> {
  // The dimension of the manifold.
  enum {
    dimension = 12
  };

  typedef kindr::minimal::HermiteTransformation<double> type;
  typedef type::QuatTransformation QuatTransformation;
  typedef type::Vector6 Vector6;
  // the "vector" typedef is used by gtsam.
  typedef Eigen::Matrix<double, dimension, 1> vector;

  // Print the type.
  static void Print(const type& T,
                    const std::string& str) {
    if(str.size() > 0) { std::cout << str << ":\n";}
    std::cout << T.getTransformation().getTransformationMatrix() << std::endl;
    std::cout << T.getTransformationDerivative() << std::endl;
  }

  // Check the equality of two values.
  static bool Equals(const type& T1,
                     const type& T2, double tol) {
    return (T1.getTransformation().getTransformationMatrix() - T2.getTransformation().getTransformationMatrix()).array().abs().maxCoeff() < tol &&
        gtsam::traits<Vector6>::Equals(T1.getTransformationDerivative(), T2.getTransformationDerivative());
  }

  // todo verify this composition of rval
  static vector Local(const type& origin, const type& other) {
    vector rval;
    rval << (other.getTransformation() * origin.getTransformation().inverted()).log(),
        (other.getTransformationDerivative() - origin.getTransformationDerivative());
    return rval;
  }

  static type Retract(const type& origin, const vector& d) {
    type rval;
    rval.setTransformation(QuatTransformation(Vector6(d.head<6>())) * origin.getTransformation());
    rval.setTransformationDerivative(Vector6(d.tail<6>()) + origin.getTransformationDerivative());
    return rval;
  }
  static int GetDimension(const type& /* origin */) {
    return dimension;
  }
};  // traits
} // namespace gtsam



namespace curves {

/// Implements the Cubic Hermite curve class. See KimKimShin paper.
/// The Hermite interpolation function is defined, with the respective Jacobians regarding  A and B:
//
/// Translations:
/// Equations for the unit interval:
// Let t_W_A, t_W_B denote the control point values (=translations) and W_v_W_A, W_v_W_B
// the control derivatives (=velocities in R³) at the interval's boundaries.
// Then (b == beta):
// p0 = t_W_A, p1 = t_W_B, p2 = W_v_W_A, p3 = W_v_W_B
// b0 = 2t³-3t²+1, b_1 = -2t³+3t², b_2 = t³-2t²+t, b_3 = t³-t²
// Spline equation:
// p(t) = p0 * b0 + p1 * b1 + p2 * b2 + p3 + b3
//
/// Rotations:
/// Equations for the unit interval:
// Let quat_W_A, quat_W_B denote the control point values (=unit quaternions) and va, vb
// the control derivatives (=angular speeds in R³) at the interval's boundaries.
// Then (w == omega, b == beta):
// w_1 = va / 3
// w_2 = log[ exp(w_1)^{-1} * quat_W_A^{-1} * quat_W_B * exp(w_3)^{-1} ]
// w_3 = vb / 3
// b_1 = t³-3t²+3t, b_2 = -2t³+3t², b_3 = t³
// Spline equation:
// q(t) = p_1 * exp(w_1*b_1) * exp(w_2*b_2) * exp(w_3*b_3)

class CubicHermiteSE3Curve : public SE3Curve {
 public:
  typedef SE3Curve::ValueType ValueType;
  typedef SE3Curve::DerivativeType DerivativeType;
  typedef kindr::minimal::HermiteTransformation<double> Coefficient;
  typedef LocalSupport2CoefficientManager<Coefficient>::TimeToKeyCoefficientMap TimeToKeyCoefficientMap;
  typedef LocalSupport2CoefficientManager<Coefficient>::CoefficientIter CoefficientIter;

  CubicHermiteSE3Curve();
  virtual ~CubicHermiteSE3Curve();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// The first valid time for the curve.
  virtual Time getMinTime() const;

  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const;

  bool isEmpty() const;

  // return number of coefficients curve is composed of
  int size() const;

  /// Extend the curve so that it can be evaluated at these times.
  /// Try to make the curve fit to the values.
  /// Note: Assumes that extend times strictly increase the curve time
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys = NULL);

  /// \brief Fit a new curve to these data points.
  ///
  /// The existing curve will be cleared.
  /// Underneath the curve should have some default policy for fitting.
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL);

  /// Evaluate the ambient space of the curve.
  virtual ValueType evaluate(Time time) const;

  /// Evaluate the curve derivatives.
  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) const;

  /// \brief Get an evaluator at this time
  virtual gtsam::Expression<ValueType> getValueExpression(const Time& time) const;

  virtual gtsam::Expression<DerivativeType> getDerivativeExpression(const Time& time, unsigned derivativeOrder) const;

  virtual void setTimeRange(Time minTime, Time maxTime);

  /// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularVelocityA(Time time);

  /// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularVelocityB(Time time);

  /// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearVelocityA(Time time);

  /// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearVelocityB(Time time);

  /// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistA(Time time);

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistB(Time time);

  /// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularDerivativeB(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearDerivativeB(unsigned derivativeOrder, Time time);

  /// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeB(unsigned derivativeOrder, Time time);

  /// Initialize a GTSAM values structure with the desired keys
  virtual void initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) const;

  /// Initialize a GTSAM values structure for all keys
  virtual void initializeGTSAMValues(gtsam::Values* values) const;

  // updates the relevant curve coefficients from the GTSAM values structure
  virtual void updateFromGTSAMValues(const gtsam::Values& values);

 private:
  LocalSupport2CoefficientManager<Coefficient> manager_;
};

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

SE3 transformationPower(SE3  T, double alpha);

SE3 composeTransformations(SE3 A, SE3 B);

SE3 inverseTransformation(SE3 T);

} // namespace curves

namespace kindr {
namespace minimal {

typedef gtsam::Expression<kindr::minimal::HermiteTransformation<double>> EHermiteTransformation;
// Expressions for HermiteTransformation -> Transformation
QuatTransformation transformationFromHermiteTransformationImplementation(
    const HermiteTransformation<double>& T, gtsam::OptionalJacobian<6, 12> HT) {
  if(HT) {
    HT->leftCols<6>().setIdentity();
    HT->rightCols<6>().setZero();
  }
  return T.getTransformation();
}

ETransformation transformationFromHermiteTransformation(
    const EHermiteTransformation& T) {
  return ETransformation(
      &transformationFromHermiteTransformationImplementation, T);
}

// Expressions for HermiteTransformation -> AngularVelocities
Eigen::Vector3d angularVelocitiesFromHermiteTransformationImplementation(
    const HermiteTransformation<double>& T, gtsam::OptionalJacobian<3, 12> HT) {
  if(HT) {
    HT->block(0,0,3,3).setZero();
    HT->block(0,3,3,3).setZero();
    HT->block(0,6,3,3).setZero();
    HT->block(0,9,3,3).setIdentity();
  }
  return T.getTransformationDerivative().tail<3>();
}

EVector3 angularVelocitiesFromHermiteTransformation(
    const EHermiteTransformation& T) {
  return EVector3(
      &angularVelocitiesFromHermiteTransformationImplementation, T);
}

// Expressions for HermiteTransformation -> Velocities
Eigen::Vector3d velocitiesFromHermiteTransformationImplementation(
    const HermiteTransformation<double>& T, gtsam::OptionalJacobian<3, 12> HT) {
  if(HT) {
    HT->block(0,0,3,3).setZero();
    HT->block(0,3,3,3).setZero();
    HT->block(0,6,3,3).setIdentity();
    HT->block(0,9,3,3).setZero();
  }
  return T.getTransformationDerivative().head<3>();
}

EVector3 velocitiesFromHermiteTransformation(
    const EHermiteTransformation& T) {
  return EVector3(
      &velocitiesFromHermiteTransformationImplementation, T);
}

/// \brief Rotate a point.
///
/// This is syntatic sugar to be able to write
/// Expression<Eigen::Vector3d> walpha = w * alpha;
/// instead of
/// Expression<Eigen::Vector3d> walpha = Expression<Eigen::Vector3d>(&vectorScaling, w, alpha);
EVector3 operator*(const EVector3& w, const double& alpha) {
  return vectorScaling(w, alpha);
}

}
}

#endif /* CURVES_CUBIC_HERMITE_SE3_CURVE_HPP */
