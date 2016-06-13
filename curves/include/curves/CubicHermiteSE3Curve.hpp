/*
 * @file CubicHermiteSE3Curve.hpp
 * @date Feb 10, 2015
 * @author Abel Gawel, Renaud Dube
 */

#ifndef CURVES_CUBIC_HERMITE_SE3_CURVE_HPP
#define CURVES_CUBIC_HERMITE_SE3_CURVE_HPP

#include "SE3Curve.hpp"
#include "LocalSupport2CoefficientManager.hpp"
#include "kindr/minimal/cubic-hermite-interpolation-gtsam.h"
#include "kindr/minimal/cubic-hermite-quaternion-gtsam.h"
#include "kindr/minimal/rotation-quaternion-gtsam.h"
#include "kindr/minimal/common-gtsam.h"
#include "gtsam/nonlinear/NonlinearFactorGraph.h"
#include "SamplingPolicy.hpp"
#include "SE3CompositionCurve.hpp"

// wrapper class for Hermite-style coefficients (made of QuatTransformation and Vector6)
namespace kindr {
namespace minimal {

typedef gtsam::Expression<kindr::minimal::RotationQuaternion> EQuaternion;
typedef gtsam::Expression<Eigen::Vector3d> EVector3;
typedef gtsam::Expression<Eigen::Matrix<double, 6, 1>> EVector6;
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
template<> struct traits<kindr::minimal::HermiteTransformation<double> > {
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
    std::cout << T.getTransformationDerivative().transpose() << std::endl;
  }

  // Check the equality of two values.
  static bool Equals(const type& T1,
                     const type& T2, double tol) {
    return (T1.getTransformation().getTransformationMatrix() - T2.getTransformation().getTransformationMatrix()).array().abs().maxCoeff() < tol &&
        gtsam::traits<Vector6>::Equals(T1.getTransformationDerivative(), T2.getTransformationDerivative());
  }

  // v_local = [log(T_w_other * T_w_origin^-1); V_other - V_origin]
  static vector Local(const type& origin, const type& other) {
    vector rval;
    rval << (other.getTransformation() * origin.getTransformation().inverted()).log(),
        (other.getTransformationDerivative() - origin.getTransformationDerivative());
    return rval;
  }

  static type Retract(const type& origin, const vector& d) {
    type rval;
    // The QuatTransformation constructor is using the exp map compatible with log
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

typedef SE3Curve::ValueType ValueType;
typedef SE3Curve::DerivativeType DerivativeType;
typedef kindr::minimal::HermiteTransformation<double> Coefficient;
typedef LocalSupport2CoefficientManager<Coefficient>::TimeToKeyCoefficientMap TimeToKeyCoefficientMap;
typedef LocalSupport2CoefficientManager<Coefficient>::CoefficientIter CoefficientIter;

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
  friend class SamplingPolicy;
 public:

  typedef kindr::minimal::HermiteTransformation<double> Coefficient;
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

  /// \brief calculate the slope between 2 coefficients
  DerivativeType calculateSlope(const Time& timeA,
                                const Time& timeB,
                                const ValueType& coeffA,
                                const ValueType& coeffB) const;

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

  void addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const;

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
  virtual void initializeGTSAMValues(gtsam::KeySet keys, gtsam::Values* values) const;

  /// Initialize a GTSAM values structure for all keys
  virtual void initializeGTSAMValues(gtsam::Values* values) const;

  // updates the relevant curve coefficients from the GTSAM values structure
  virtual void updateFromGTSAMValues(const gtsam::Values& values);

  // set the minimum sampling period
  void setMinSamplingPeriod(Time time);

  /// \brief Set the sampling ratio.
  ///   eg. 4 will add a coefficient every 4 extend
  void setSamplingRatio(const int ratio);

  // clear the curve
  virtual void clear();

  /// \brief Perform a rigid transformation on the left side of the curve
  void transformCurve(const ValueType T);

  virtual Time getTimeAtKey(gtsam::Key key) const;

  void saveCurveTimesAndValues(const std::string& filename) const;

  void saveCurveAtTimes(const std::string& filename, std::vector<Time> times) const;

  void saveCorrectionCurveAtTimes(const std::string& filename, std::vector<Time> times) const {};

  void getCurveTimes(std::vector<Time>* outTimes) const;

  /// \brief Returns the number of coefficients in the correction curve
  int correctionSize() const {return 0;};

  /// \brief Fold in the correction curve into the base curve and reinitialize
  ///        correction curve coefficients to identity transformations.
  void foldInCorrections() {};

  /// \brief Add coefficients to the correction curve at given times.
  void setCorrectionTimes(const std::vector<Time>& times) {};

  /// \brief Remove a correction coefficient at the specified time.
  void removeCorrectionCoefficientAtTime(Time time) {};

  /// \brief Set the correction coefficient value at the specified time.
  void setCorrectionCoefficientAtTime(Time time, ValueType value) {};

  /// \brief Reset the correction curve to identity values with knots at desired times
  void resetCorrectionCurve(const std::vector<Time>& times) {};

  /// \brief Set the base curve to given values with knots at desired times
  /// Resets the curve beforehand.
  void setBaseCurve(const std::vector<Time>& times, const std::vector<ValueType>& values) {};

  /// \brief Add / replace the given coefficients without resetting the curve.
  void setBaseCurvePart(const std::vector<Time>& times, const std::vector<ValueType>& values) {};

  /// \brief Modifies values of the base coefficient in batch, starting at times[0] and assuming that
  /// a coefficient exists at all the specified times.
  void modifyBaseCoefficientsValuesInBatch(const std::vector<Time>& times, const std::vector<ValueType>& values) {};

  void getBaseCurveTimes(std::vector<Time>* outTimes) const {};

  void getBaseCurveTimesInWindow(std::vector<Time>* outTimes, Time begTime, Time endTime) const {};

  // return number of coefficients curve is composed of
  int baseSize() const {return size();};

  void saveCorrectionCurveTimesAndValues(const std::string& filename) const {};
 private:
  LocalSupport2CoefficientManager<Coefficient> manager_;
  SamplingPolicy hermitePolicy_;
};

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

SE3 transformationPower(SE3  T, double alpha);

SE3 composeTransformations(SE3 A, SE3 B);

SE3 inverseTransformation(SE3 T);


// implements the special (extend) policies for Cubic Hermite curves
template <>
inline Key SamplingPolicy::defaultExtend<CubicHermiteSE3Curve, ValueType>(const Time& time,
                  const ValueType& value,
                  CubicHermiteSE3Curve* curve) {

  DerivativeType derivative;
  //cases:

  // manager is empty
  if (curve->manager_.size() == 0) {
    derivative << 0,0,0,0,0,0;
    // 1 value in manager (2 in total)
  } else if (curve->manager_.size() == 1) {
    // get latest coefficient from manager
    CoefficientIter last = curve->manager_.coefficientBegin();
    // calculate slope
    // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
    derivative = curve->calculateSlope(last->first,
                                time,
                                last->second.coefficient.getTransformation(),
                                value);
    // update previous coefficient
    Coefficient updated(last->second.coefficient.getTransformation(), derivative);
    curve->manager_.updateCoefficientByKey(last->second.key, updated);
    // more than 1 values in manager
  } else if (curve->manager_.size() > 1) {
    // get latest 2 coefficients from manager
    CoefficientIter rVal0, rVal1;
    CoefficientIter last = --curve->manager_.coefficientEnd();
    curve->manager_.getCoefficientsAt(last->first,
                               &rVal0,
                               &rVal1);

    // update derivative of previous coefficient
    DerivativeType derivative0;
    derivative0 << curve->calculateSlope(rVal0->first,
                                  time,
                                  rVal0->second.coefficient.getTransformation(),
                                  value);
    Coefficient updated(rVal1->second.coefficient.getTransformation(), derivative0);
    curve->manager_.updateCoefficientByKey(rVal1->second.key, updated);

    // calculate slope
    derivative << curve->calculateSlope(rVal1->first,
                                 time,
                                 rVal1->second.coefficient.getTransformation(),
                                 value);
  }
  measurementsSinceLastExtend_ = 0;
  lastExtend_ = time;
  return curve->manager_.insertCoefficient(time, Coefficient(value, derivative));
}

template <>
inline Key SamplingPolicy::interpolationExtend<CubicHermiteSE3Curve, ValueType>(const Time& time,
                        const ValueType& value,
                        CubicHermiteSE3Curve* curve) {
  DerivativeType derivative;
  if (measurementsSinceLastExtend_ == 0) {
    // extend curve with new interpolation coefficient if necessary
    CoefficientIter rValInterp = --curve->manager_.coefficientEnd();
    CoefficientIter last = --curve->manager_.coefficientEnd();
    derivative << last->second.coefficient.getTransformationDerivative();
  } else {
    // assumes the interpolation coefficient is already set (at end of curve)
    // assumes same velocities as last Coefficient
    CoefficientIter rVal0, rValInterp;
    CoefficientIter last = --curve->manager_.coefficientEnd();
    curve->manager_.getCoefficientsAt(last->first,
                               &rVal0,
                               &rValInterp);

    derivative << curve->calculateSlope(rVal0->first,
                                 time,
                                 rVal0->second.coefficient.getTransformation(),
                                 value);

    // update the interpolated coefficient with given values and velocities from last coefficeint
    curve->manager_.removeCoefficientAtTime(rValInterp->first);
  }

  ++measurementsSinceLastExtend_;
  return curve->manager_.insertCoefficient(time, Coefficient(value, derivative));
}

template<>
inline void SamplingPolicy::extend<CubicHermiteSE3Curve, ValueType>(const std::vector<Time>& times,
                                                             const std::vector<ValueType>& values,
                                                             CubicHermiteSE3Curve* curve,
                                                             std::vector<Key>* outKeys) {

  for (int i = 0; i < times.size(); ++i) {
    // ensure time strictly increases
    CHECK((times[i] > curve->manager_.getMaxTime()) || curve->manager_.size() == 0) << "curve can only be extended into the future. Requested = "
        << times[i] << " < curve max time = " << curve->manager_.getMaxTime();
    if (curve->manager_.size() == 0) {
      defaultExtend(times[i], values[i], curve);
    } else if((measurementsSinceLastExtend_ >= minimumMeasurements_ &&
        lastExtend_ + minSamplingPeriod_ < times[i])) {
      // delete interpolated coefficient
      CoefficientIter last = --curve->manager_.coefficientEnd();
      curve->manager_.removeCoefficientAtTime(last->first);
      // todo write outkeys
      defaultExtend(times[i], values[i], curve);
    } else {
      interpolationExtend(times[i], values[i], curve);
    }
  }
}

} // namespace curves

namespace kindr {
namespace minimal {

typedef gtsam::Expression<kindr::minimal::HermiteTransformation<double>> EHermiteTransformation;
typedef curves::SE3Curve::ValueType ValueType;
typedef curves::SE3Curve::DerivativeType DerivativeType;
typedef kindr::minimal::HermiteTransformation<double> Coefficient;
// Expressions for HermiteTransformation -> Transformation
static QuatTransformation transformationFromHermiteTransformationImplementation(
    const HermiteTransformation<double>& T, gtsam::OptionalJacobian<6, 12> HT) {
  if(HT) {
    HT->leftCols<6>().setIdentity();
    HT->rightCols<6>().setZero();
  }
  return T.getTransformation();
}

static ETransformation transformationFromHermiteTransformation(
    const EHermiteTransformation& T) {
  return ETransformation(
      &transformationFromHermiteTransformationImplementation, T);
}

// Expressions for HermiteTransformation -> AngularVelocities
static Eigen::Vector3d angularVelocitiesFromHermiteTransformationImplementation(
    const HermiteTransformation<double>& T, gtsam::OptionalJacobian<3, 12> HT) {
  if(HT) {
    HT->block(0,0,3,3).setZero();
    HT->block(0,3,3,3).setZero();
    HT->block(0,6,3,3).setZero();
    HT->block(0,9,3,3).setIdentity();
  }
  return T.getTransformationDerivative().tail<3>();
}

static EVector3 angularVelocitiesFromHermiteTransformation(
    const EHermiteTransformation& T) {
  return EVector3(
      &angularVelocitiesFromHermiteTransformationImplementation, T);
}

// Expressions for HermiteTransformation -> Velocities
static Eigen::Vector3d velocitiesFromHermiteTransformationImplementation(
    const HermiteTransformation<double>& T, gtsam::OptionalJacobian<3, 12> HT) {
  if(HT) {
    HT->block(0,0,3,3).setZero();
    HT->block(0,3,3,3).setZero();
    HT->block(0,6,3,3).setIdentity();
    HT->block(0,9,3,3).setZero();
  }
  return T.getTransformationDerivative().head<3>();
}

static EVector3 velocitiesFromHermiteTransformation(
    const EHermiteTransformation& T) {
  return EVector3(
      &velocitiesFromHermiteTransformationImplementation, T);
}

static Coefficient composeHermiteTransformationImplementation(
    const ValueType& T, const DerivativeType& v,
    gtsam::OptionalJacobian<12, 6> H1,
    gtsam::OptionalJacobian<12, 6> H2) {
  if (H1) {
    H1->topRows<6>().setIdentity();
    H1->bottomRows<6>().setZero();
  }
  if (H2) {
    H2->topRows<6>().setZero();
    H2->bottomRows<6>().setIdentity();
  }
  return Coefficient(T, v);
}

// Expression for Transformation & Velocities -> HermiteTransformation
static EHermiteTransformation composeHermiteTransformation(
    const ETransformation& T,
    const EVector6& v) {
  return EHermiteTransformation(
      &composeHermiteTransformationImplementation, T, v);
}
/// \brief Scale a point.
///
/// This is syntatic sugar to be able to write
/// Expression<Eigen::Vector3d> walpha = w * alpha;
/// instead of
/// Expression<Eigen::Vector3d> walpha = Expression<Eigen::Vector3d>(&vectorScaling, w, alpha);
static EVector3 operator*(const EVector3& w, const double& alpha) {
  return vectorScaling(w, alpha);
}

}
}


#endif /* CURVES_CUBIC_HERMITE_SE3_CURVE_HPP */
