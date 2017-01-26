/*
 * CubicHermiteSE3Curve.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: Abel Gawel, Renaud Dube, Péter Fankhauser, Christian Gehring
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <iostream>

#include "curves/CubicHermiteSE3Curve.hpp"
#include "curves/SlerpSE3Curve.hpp"

namespace curves {

CubicHermiteSE3Curve::CubicHermiteSE3Curve() : SE3Curve() {
  hermitePolicy_.setMinimumMeasurements(4);
}

CubicHermiteSE3Curve::~CubicHermiteSE3Curve() {}

void CubicHermiteSE3Curve::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "======= Cubic Hermite SE3 CURVE =========" << std::endl;
  std::cout << str << std::endl;
  std::cout << "num of coefficients: " << manager_.size() << std::endl;
  std::cout << "dimension: " << 6 << std::endl;
  std::stringstream ss;
  std::vector<Key> keys;
  std::vector<Time> times;
  manager_.getTimes(&times);
  manager_.getKeys(&keys);
  std::cout << "curve defined between times: " << manager_.getMinTime() <<
      " and " << manager_.getMaxTime() <<std::endl;
  std::cout <<"=========================================" <<std::endl;
  for (size_t i = 0; i < manager_.size(); i++) {
    ss << "coefficient " << keys[i] << ": ";
    ss << manager_.getCoefficientByKey(keys[i]).getTransformation() << std::endl;
    ss << " | time: " << times[i];
    ss << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
}

bool CubicHermiteSE3Curve::writeEvalToFile(const std::string& filename, int nSamples) const {
  FILE* fp = fopen(filename.c_str(), "w");
  if (fp==NULL) {
    std::cout << "Could not open file to write" << std::endl;
    return false;
  }
  fprintf(fp, "t ");
  fprintf(fp, "px py pz ");
  fprintf(fp, "rw rx ry rz ");
  fprintf(fp, "vx vy vz ");
  fprintf(fp, "wx wy wz ");
  fprintf(fp, "\n");

  Time dt = (getMaxTime()-getMinTime())/(nSamples-1);
  ValueType pose;
  DerivativeType twist;
  for (Time t = getMinTime(); t < getMaxTime(); t+=dt) {
    if(!evaluate(pose, t)) {
      std::cout << "Could not evaluate at time " << t << std::endl;
      fclose(fp);
      return false;
    }
    if(!evaluateDerivative(twist, t, 1)) {
      std::cout << "Could not evaluate derivative at time " << t << std::endl;
      fclose(fp);
      return false;
    }
    fprintf(fp, "%lf ", t);
    fprintf(fp, "%lf %lf %lf ", pose.getPosition().x(), pose.getPosition().y(), pose.getPosition().z());
    fprintf(fp, "%lf %lf %lf %lf ", pose.getRotation().w(), pose.getRotation().x(), pose.getRotation().y(), pose.getRotation().z());
    fprintf(fp, "%lf %lf %lf ", twist.getTranslationalVelocity().x(), twist.getTranslationalVelocity().y(), twist.getTranslationalVelocity().z());
    fprintf(fp, "%lf %lf %lf ", twist.getRotationalVelocity().x(), twist.getRotationalVelocity().y(), twist.getRotationalVelocity().z());
    fprintf(fp, "\n");
  }
  fclose(fp);

  return true;
}

Time CubicHermiteSE3Curve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time CubicHermiteSE3Curve::getMinTime() const {
  return manager_.getMinTime();
}

bool CubicHermiteSE3Curve::isEmpty() const {
  std::vector<Time> outTimes;
  manager_.getTimes(&outTimes);
  return outTimes.empty();
}

int CubicHermiteSE3Curve::size() const {
  return manager_.size();
}

void CubicHermiteSE3Curve::fitCurve(const std::vector<Time>& times,
                                    const std::vector<ValueType>& values,
                                    std::vector<Key>* outKeys) {
  fitCurveWithDerivatives(times, values, DerivativeType(), DerivativeType(), outKeys);
}

void CubicHermiteSE3Curve::fitCurveWithDerivatives(const std::vector<Time>& times,
                                    const std::vector<ValueType>& values,
                                    const DerivativeType& initialDerivative,
                                    const DerivativeType& finalDerivative,
                                    std::vector<Key>* outKeys)
{
  CHECK_EQ(times.size(), values.size());
  clear();

  // construct the Hemrite coefficients
  std::vector<Coefficient> coefficients;
  // fill the coefficients with ValueType and DerivativeType
  // use Catmull-Rom interpolation for derivatives on knot points
  for (size_t i = 0; i < times.size(); ++i) {
    DerivativeType derivative;
    // catch the boundaries (i == 0 && i == max)
    if (i == 0) {
      // First key.
      if (times.size() > 1) {
//        derivative = calculateSlope(times[0], times[1], values[0], values[1]);
        derivative = initialDerivative;
      } else {
        // set velocities == 0 for start point if only one coefficient
        derivative = initialDerivative;
      }
    } else if (i == times.size() - 1) {
      // Last key.
//      derivative = calculateSlope(times[i-1], times[i], values[i-1], values[i]);
      derivative = finalDerivative;
    } else {
      // Other keys.
      derivative = calculateSlope(times[i-1], times[i+1], values[i-1], values[i+1]);
    }
    coefficients.push_back(Coefficient(values[i], derivative));
  }

  manager_.insertCoefficients(times, coefficients, outKeys);
}

void CubicHermiteSE3Curve::fitPeriodicCurve(const std::vector<Time>& times,
                                           const std::vector<ValueType>& values,
                                           std::vector<Key>* outKeys)
{
  /* We assume that the first and last points coincide.
   *
   */
  const size_t nPoints = times.size();
  CubicHermiteSE3Curve::DerivativeType derivative = calculateSlope(times[nPoints-2], times[1], values[nPoints-2], values[1]);
  fitCurveWithDerivatives(times, values, derivative, derivative, outKeys);
}

CubicHermiteSE3Curve::DerivativeType CubicHermiteSE3Curve::calculateSlope(const Time& timeA,
                                                                          const Time& timeB,
                                                                          const ValueType& T_W_A,
                                                                          const ValueType& T_W_B) const {
  const double inverse_dt_sec = 1.0/double(timeB - timeA);
  // Original curves implementation was buggy for 180 deg flips.

  // Calculate the global angular velocity:
  const Eigen::Vector3d angularVelocity_rad_s = T_W_B.getRotation().boxMinus(T_W_A.getRotation()) * inverse_dt_sec;
  const Eigen::Vector3d velocity_m_s = (T_W_B.getPosition().vector() - T_W_A.getPosition().vector()) * inverse_dt_sec;
  // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
  return DerivativeType(velocity_m_s, angularVelocity_rad_s);
}

void CubicHermiteSE3Curve::extend(const std::vector<Time>& times,
                                  const std::vector<ValueType>& values,
                                  std::vector<Key>* outKeys) {

  // New values in extend first need to be checked if they can be added to curve
  // otherwise the most recent coefficient will be an interpolation based on the last
  // coefficient and the requested time-value pair

  // switch extension curve between default and interpolation (with regard to policy)
  // - default extend if the minimum time gap between 2 coefficients is satisfied
  // - default extend if enough measurements were taken between 2 coefficients, so
  //   the curve is well defined
  // - default extend if curve is empty
  // - interpolation extend otherwise

//  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
//  hermitePolicy_.extend<CubicHermiteSE3Curve, ValueType>(times, values, this, outKeys);
  throw std::runtime_error("CubicHermiteSE3Curve::extend is not implemented!");
}


bool CubicHermiteSE3Curve::evaluate(ValueType& value, Time time) const {
  // Check if the curve is only defined at this one time
  if (manager_.getMaxTime() == time && manager_.getMinTime() == time) {
    value =  manager_.coefficientBegin()->second.coefficient.getTransformation();
    return true;
  }
  else {
    CoefficientIter a, b;
    bool success = manager_.getCoefficientsAt(time, &a, &b);
    if(!success) {
      std::cerr << "Unable to get the coefficients at time " << time << std::endl;
      return false;
    }

    // read out transformation from coefficient
    const SE3 T_W_A = a->second.coefficient.getTransformation();
    const SE3 T_W_B = b->second.coefficient.getTransformation();

    // read out derivative from coefficient
    const Twist d_W_A = a->second.coefficient.getTransformationDerivative();
    const Twist d_W_B = b->second.coefficient.getTransformationDerivative();

    // make alpha
    const double dt_sec = (b->first - a->first);// * 1e-9;
    const double alpha = double(time - a->first)/(b->first - a->first);

    // Implemantation of Hermite Interpolation not easy and not fun (without expressions)!

    // translational part (easy):
    const double alpha2 = alpha * alpha;
    const double alpha3 = alpha2 * alpha;

    const double beta0 = 2.0 * alpha3 - 3.0 * alpha2 + 1.0;
    const double beta1 = -2.0 * alpha3 + 3.0 * alpha2;
    const double beta2 = alpha3 - 2.0 * alpha2 + alpha;
    const double beta3 = alpha3 - alpha2;

    /**************************************************************************************
     *  Translational part:
     **************************************************************************************/
    const SE3::Position translation(T_W_A.getPosition().vector() * beta0
                                  + T_W_B.getPosition().vector() * beta1
                                  + d_W_A.getTranslationalVelocity().vector() * (beta2 * dt_sec)
                                  + d_W_B.getTranslationalVelocity().vector() * (beta3 * dt_sec));

    /**************************************************************************************
     *  Rotational part:
     **************************************************************************************/
    const double dt_sec_third = dt_sec / 3.0;
    const Eigen::Vector3d scaled_d_W_A = dt_sec_third * d_W_A.getRotationalVelocity().vector();
    const Eigen::Vector3d scaled_d_W_B = dt_sec_third * d_W_B.getRotationalVelocity().vector();

    // d_W_A contains the global angular velocity, but we need the local angular velocity.
    const Eigen::Vector3d w1 = T_W_A.getRotation().inverseRotate(scaled_d_W_A);
    const Eigen::Vector3d w3 = T_W_B.getRotation().inverseRotate(scaled_d_W_B);
    const RotationQuaternion expW1_inv = RotationQuaternion().exponentialMap(-w1);
    const RotationQuaternion expW3_inv = RotationQuaternion().exponentialMap(-w3);
    const RotationQuaternion expW1_Inv_qWB_expW3 = expW1_inv * T_W_A.getRotation().inverted() * T_W_B.getRotation() * expW3_inv;
    const Eigen::Vector3d w2 = expW1_Inv_qWB_expW3.logarithmicMap();

    const double dBeta1 = alpha3 - 3.0 * alpha2 + 3.0 * alpha;
    const double dBeta2 = -2.0 * alpha3 + 3.0 * alpha2;
    const double dBeta3 = alpha3;

    const SO3 w1_dBeta1_exp = RotationQuaternion().exponentialMap(dBeta1 * w1);
    const SO3 w2_dBeta2_exp = RotationQuaternion().exponentialMap(dBeta2 * w2);
    const SO3 w3_dBeta3_exp = RotationQuaternion().exponentialMap(dBeta3 * w3);

    const RotationQuaternion rotation = T_W_A.getRotation() * w1_dBeta1_exp * w2_dBeta2_exp * w3_dBeta3_exp;

    value = SE3(translation, rotation);
    return true;
  }
  return false;
}

bool CubicHermiteSE3Curve::evaluateDerivative(DerivativeType& derivative,
    Time time, unsigned int derivativeOrder) const
{
  if (derivativeOrder == 1) {
    // Check if the curve is only defined at this one time
    if (manager_.getMaxTime() == time && manager_.getMinTime() == time) {
      derivative = manager_.coefficientBegin()->second.coefficient.getTransformationDerivative();
      return true;
    }
    else {
      CoefficientIter a, b;
      bool success = manager_.getCoefficientsAt(time, &a, &b);
      if(!success) {
        std::cerr << "Unable to get the coefficients at time " << time << std::endl;
        return false;
      }

      // read out transformation from coefficient
      const SE3 T_W_A = a->second.coefficient.getTransformation();
      const SE3 T_W_B = b->second.coefficient.getTransformation();

      // read out derivative from coefficient
      const Twist d_W_A = a->second.coefficient.getTransformationDerivative();
      const Twist d_W_B = b->second.coefficient.getTransformationDerivative();

      // make alpha
      double dt_sec = (b->first - a->first);
      const double one_over_dt_sec = 1.0/dt_sec;
      double alpha = double(time - a->first)/dt_sec;

      const double alpha2 = alpha * alpha;
      const double alpha3 = alpha2 * alpha;

      /**************************************************************************************
       *  Translational part:
       **************************************************************************************/
      // Implementation of translation
      const double gamma0 = 6.0*(alpha2 - alpha);
      const double gamma1 = 3.0*alpha2 - 4.0*alpha + 1.0;
      const double gamma2 = 6.0*(alpha - alpha2);
      const double gamma3 = 3.0*alpha2 - 2.0*alpha;

      const Eigen::Vector3d velocity_m_s = T_W_A.getPosition().vector()*(gamma0*one_over_dt_sec)
                                         + d_W_A.getTranslationalVelocity().vector()*(gamma1)
                                         + T_W_B.getPosition().vector()*(gamma2*one_over_dt_sec)
                                         + d_W_B.getTranslationalVelocity().vector()*(gamma3);


      /**************************************************************************************
       *  Rotational part:
       **************************************************************************************/
      const double one_minus_alpha = (1.0 - alpha);
      const double one_minus_alpha_2 = one_minus_alpha * one_minus_alpha;
      const double one_minus_alpha_3 = one_minus_alpha * one_minus_alpha_2;

      const double beta1 = 1.0 - one_minus_alpha_3;
      const double dbeta1 = 3.0*one_minus_alpha_2;
      const double beta2 = 3.0*alpha2 - 2.0*alpha3;
      const double dbeta2 = 6.0*alpha*one_minus_alpha;
      const double beta3 = alpha3;
      const double dbeta3 = 3.0*alpha2;

      const double one_third = 1.0 / 3.0;
      const Eigen::Vector3d scaled_d_W_A = (one_third*dt_sec ) * d_W_A.getRotationalVelocity().vector();
      const Eigen::Vector3d scaled_d_W_B = (one_third*dt_sec ) * d_W_B.getRotationalVelocity().vector();

      const Eigen::Vector3d w1 = T_W_A.getRotation().inverseRotate(scaled_d_W_A);
      const Eigen::Vector3d w3 = T_W_B.getRotation().inverseRotate(scaled_d_W_B);
      const RotationQuaternion expW1_inv = RotationQuaternion().exponentialMap(-w1);
      const RotationQuaternion expW3_inv = RotationQuaternion().exponentialMap(-w3);

      const RotationQuaternion expW1_Inv_qWB_expW3 = expW1_inv * T_W_A.getRotation().inverted() * T_W_B.getRotation() * expW3_inv;

      const Eigen::Vector3d w2 = expW1_Inv_qWB_expW3.logarithmicMap();

      const SO3 w1_beta1_exp = RotationQuaternion().exponentialMap((beta1) * w1);
      const SO3 w2_beta2_exp = RotationQuaternion().exponentialMap((beta2) * w2);
      const SO3 w3_beta3_exp = RotationQuaternion().exponentialMap((beta3) * w3);

      const RotationQuaternion w1_dbeta1(0.0, dbeta1 * w1);
      const RotationQuaternion w2_dbeta2(0.0, dbeta2 * w2);
      const RotationQuaternion w3_dbeta3(0.0, dbeta3 * w3);

      const Eigen::Vector4d diff =    ((T_W_A.getRotation() * w1_beta1_exp * w1_dbeta1    * w2_beta2_exp * w3_beta3_exp).vector()
                              + (T_W_A.getRotation() * w1_beta1_exp * w2_beta2_exp * w2_dbeta2    * w3_beta3_exp).vector()
                              + (T_W_A.getRotation() * w1_beta1_exp * w2_beta2_exp * w3_beta3_exp * w3_dbeta3   ).vector())*one_over_dt_sec;

      const RotationQuaternion qDiff(diff);
      ValueType q;
      if(!evaluate(q, time)) {
        return false;
      }
      // This is the global angular velocity
      const Eigen::Vector3d angularVelocity_rad_s = q.getRotation().rotate((q.getRotation().inverted()*qDiff).imaginary());

      // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries

      derivative = DerivativeType(velocity_m_s, angularVelocity_rad_s);
      return true;
    }
  }
  else {
    std::cerr << "CubicHermiteSE3Curve::evaluateDerivative: higher order derivatives are not implemented!";
    return false;
  }
}

bool CubicHermiteSE3Curve::evaluateLinearAcceleration(kindr::Acceleration3D& linearAcceleration, Time time) {

  CoefficientIter a, b;
  bool success = manager_.getCoefficientsAt(time, &a, &b);
  if(!success) {
    std::cerr << "Unable to get the coefficients at time " << time << std::endl;
    return false;
  }

  // read out transformation from coefficient
  const SE3 T_W_A = a->second.coefficient.getTransformation();
  const SE3 T_W_B = b->second.coefficient.getTransformation();

  // read out derivative from coefficient
  const Twist d_W_A = a->second.coefficient.getTransformationDerivative();
  const Twist d_W_B = b->second.coefficient.getTransformationDerivative();

  // make alpha
  double dt_sec = (b->first - a->first);
  const double one_over_dt_sec = 1.0/dt_sec;
  double alpha = double(time - a->first)/dt_sec;
  const double d_alpha = one_over_dt_sec;

  /**************************************************************************************
   *  Translational part:
   **************************************************************************************/
  // Implementation of translation
  const double d_gamma0 = 6.0*(2*alpha - 1.0)*d_alpha;
  const double d_gamma1 = (6.0*alpha - 4.0)*d_alpha;
  const double d_gamma2 = 6.0*(1.0 - 2.0*alpha)*d_alpha;
  const double d_gamma3 = (6.0*alpha - 2.0)*d_alpha;

  linearAcceleration = kindr::Acceleration3D(T_W_A.getPosition().vector()*d_gamma0*one_over_dt_sec + d_W_A.getTranslationalVelocity().vector()*d_gamma1 +
                                             T_W_B.getPosition().vector()*d_gamma2*one_over_dt_sec + d_W_B.getTranslationalVelocity().vector()*d_gamma3);


  return true;
}

/// \brief forms cubic Hermite interpolation into a binary expression with 2 leafs and binds alpha into it,
///        uses break down of expression into its operations
//gtsam::Expression<typename CubicHermiteSE3Curve::ValueType>
//CubicHermiteSE3Curve::getValueExpression(const Time& time) const {
//
//  using namespace gtsam;
//
//  // CoefficientType is HermiteTransformation
//  CoefficientIter rval0, rval1;
//  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
//  CHECK(success) << "Unable to get the coefficients at time " << time;
//
//  // shall return the Transformation Expression for Hermite
//  EHermiteTransformation leaf1(rval0->second.key);
//  EHermiteTransformation leaf2(rval1->second.key);
//
//  double dt_sec = (rval1->first - rval0->first) * 1e-9;
//  double alpha = double(time - rval0->first)/(rval1->first - rval0->first);
//
//  // construct the necessary Expressions (and underlying Jacobians) for evaluation
//  ETransformation transformation1 = kindr::minimal::transformationFromHermiteTransformation(leaf1);
//  EVector3 angVel1 = kindr::minimal::angularVelocitiesFromHermiteTransformation(leaf1);
//  EVector3 vel1 = kindr::minimal::velocitiesFromHermiteTransformation(leaf1);
//  EQuaternion quat1 = kindr::minimal::rotationFromTransformation(transformation1);
//  EVector3 trans1 = kindr::minimal::translationFromTransformation(transformation1);
//
//  ETransformation transformation2 = kindr::minimal::transformationFromHermiteTransformation(leaf2);
//  EVector3 angVel2 = kindr::minimal::angularVelocitiesFromHermiteTransformation(leaf2);
//  EVector3 vel2 = kindr::minimal::velocitiesFromHermiteTransformation(leaf2);
//  EQuaternion quat2 = kindr::minimal::rotationFromTransformation(transformation2);
//  EVector3 trans2 = kindr::minimal::translationFromTransformation(transformation2);
//
//  EQuaternion quat = kindr::minimal::hermiteQuaternionInterpolation(quat1,
//                                                                    angVel1,
//                                                                    quat2,
//                                                                    angVel2,
//                                                                    alpha,
//                                                                    dt_sec);
//
//  EVector3 trans = kindr::minimal::hermiteInterpolation<Vector3>(trans1,
//                                                                 vel1,
//                                                                 trans2,
//                                                                 vel2,
//                                                                 alpha,
//                                                                 dt_sec);
//
//  return kindr::minimal::transformationFromComponents(quat, trans);
//}

//gtsam::Expression<typename CubicHermiteSE3Curve::DerivativeType>
//CubicHermiteSE3Curve::getDerivativeExpression(const Time& time, unsigned derivativeOrder) const {
//  // \todo Abel and Renaud
//  CHECK(false) << "Not implemented";
//}


//void CubicHermiteSE3Curve::addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const {
//
//  gtsam::noiseModel::Constrained::shared_ptr priorNoise = gtsam::noiseModel::Constrained::All(gtsam::traits<Coefficient>::dimension);
//
//  //Add one fixed prior at priorTime and two before to ensure that at least two
//  CoefficientIter rVal0, rVal1;
//  manager_.getCoefficientsAt(priorTime, &rVal0, &rVal1);
//
//  gtsam::ExpressionFactor<Coefficient> f0(priorNoise,
//                                          rVal0->second.coefficient,
//                                          gtsam::Expression<Coefficient>(rVal0->second.key));
//  gtsam::ExpressionFactor<Coefficient> f1(priorNoise,
//                                          rVal1->second.coefficient,
//                                          gtsam::Expression<Coefficient>(rVal1->second.key));
//  graph->push_back(f0);
//  graph->push_back(f1);
//
//}

void CubicHermiteSE3Curve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

/// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateTwistA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateTwistB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}

void CubicHermiteSE3Curve::setMinSamplingPeriod(Time time) {
  hermitePolicy_.setMinSamplingPeriod(time);
}

///   eg. 4 will add a coefficient every 4 extend
void CubicHermiteSE3Curve::setSamplingRatio(const int ratio) {
  hermitePolicy_.setMinimumMeasurements(ratio);
}

void CubicHermiteSE3Curve::clear() {
  manager_.clear();
}

void CubicHermiteSE3Curve::transformCurve(const ValueType T) {
  //todo
}

void CubicHermiteSE3Curve::saveCurveTimesAndValues(const std::string& filename) const {
  std::vector<Time> curveTimes;
  manager_.getTimes(&curveTimes);

  saveCurveAtTimes(filename, curveTimes);
}

void CubicHermiteSE3Curve::saveCurveAtTimes(const std::string& filename, std::vector<Time> times) const {
  Eigen::VectorXd v(7);

  std::vector<Eigen::VectorXd> curveValues;
  ValueType val;
  for (size_t i = 0; i < times.size(); ++i) {
    evaluate(val, times[i]);
    v << val.getPosition().x(), val.getPosition().y(), val.getPosition().z(),
        val.getRotation().w(), val.getRotation().x(), val.getRotation().y(), val.getRotation().z();
    curveValues.push_back(v);
  }

  writeTimeVectorCSV(filename, times, curveValues);
}

void CubicHermiteSE3Curve::getCurveTimes(std::vector<Time>* outTimes) const {
  manager_.getTimes(outTimes);
}

} // namespace curves
