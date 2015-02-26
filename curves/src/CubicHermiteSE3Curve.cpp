/*
 * @file CubicHermiteSE3Curve.cpp
 * @date Feb 10, 2015
 * @author Abel Gawel, Renaud Dube
 */

#include <curves/CubicHermiteSE3Curve.hpp>
#include "gtsam/nonlinear/ExpressionFactor.h"
#include <iostream>

namespace curves {

CubicHermiteSE3Curve::CubicHermiteSE3Curve() : SE3Curve() {}

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
    gtsam::traits<Coefficient>::Print(manager_.getCoefficientByKey(keys[i]),ss.str());
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
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

typedef gtsam::Expression<kindr::minimal::RotationQuaternion> EQuaternion;
typedef gtsam::Expression<Eigen::Vector3d> EVector3;
typedef gtsam::Expression<kindr::minimal::QuatTransformation> ETransformation;
typedef gtsam::Expression<kindr::minimal::HermiteTransformation<double>> EHermiteTransformation;


void CubicHermiteSE3Curve::fitCurve(const std::vector<Time>& times,
                                    const std::vector<ValueType>& values,
                                    std::vector<Key>* outKeys) {
  CHECK_EQ(times.size(), values.size());
  CHECK_GE(times.size(), 3) << "Hermite curve must be defined by > 2 coefficients";

  // construct the Hemrite coefficients
  std::vector<Coefficient> coefficients;
  // fill the coefficients with ValueType and DerivativeType
  // use Catmull-Rom interpolation for derivatives on knot points
  for (size_t i = 0; i < times.size(); ++i) {

    DerivativeType derivative;
    // catch the boundaries (i == 0 && i == max)
    if (i == 0) {

      if (times.size() > 1) {

        Time period = times[1] - times[0];
        SE3 T_A_B = kindr::minimal::invertAndComposeImplementation(values[0], values[1], boost::none, boost::none);
        SO3 q_A_B = T_A_B.getRotation();
        kindr::minimal::AngleAxisTemplate<double> aa_A_B(q_A_B);
        double angle = aa_A_B.angle();
        Eigen::Vector3d angularVelocity = aa_A_B.axis() * angle/period * 1e9;
        Eigen::Vector3d velocity = T_A_B.getPosition() * 1/period * 1e9;
        // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
        derivative << velocity, angularVelocity;

      } else {

        // set velocities == 0 for start point if only one coefficient
        derivative << 0,0,0,0,0,0;
      }

    } else if (i == times.size() - 1) {

      Time period = times[i] - times[i-1];
      SE3 T_A_B = kindr::minimal::invertAndComposeImplementation(values[i-1], values[i], boost::none, boost::none);
      SO3 q_A_B = T_A_B.getRotation();
      kindr::minimal::AngleAxisTemplate<double> aa_A_B(q_A_B);
      double angle = aa_A_B.angle();
      Eigen::Vector3d angularVelocity = aa_A_B.axis() * angle/period * 1e9;
      Eigen::Vector3d velocity = T_A_B.getPosition() * 1/period * 1e9;
      // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
      derivative << velocity, angularVelocity;

    } else {

      Time period = times[i+1] - times[i-1];
      SE3 T_A_B = kindr::minimal::invertAndComposeImplementation(values[i-1], values[i+1], boost::none, boost::none);
      SO3 q_A_B = T_A_B.getRotation();
      kindr::minimal::AngleAxisTemplate<double> aa_A_B(q_A_B);
      double angle = aa_A_B.angle();
      Eigen::Vector3d angularVelocity = aa_A_B.axis() * angle/period * 1e9;
      Eigen::Vector3d velocity = T_A_B.getPosition() * 1/period * 1e9;
      // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
      derivative << velocity, angularVelocity;
    }

    coefficients.push_back(Coefficient(values[i], derivative));
  }

  manager_.insertCoefficients(times, coefficients, outKeys);
}

CubicHermiteSE3Curve::DerivativeType CubicHermiteSE3Curve::calculateSlope(const Time& timeA_ns,
                                                                          const Time& timeB_ns,
                                                                          const ValueType& T_W_A,
                                                                          const ValueType& T_W_B) const {
  double inverse_dt_sec = 1e9/double(timeB_ns - timeA_ns);
  SO3 R_W_A_B = T_W_B.getRotation() * T_W_A.getRotation().inverted();
  kindr::minimal::AngleAxisTemplate<double> aa_W_A_B(R_W_A_B);
  double angle = aa_W_A_B.angle();
  Eigen::Vector3d angularVelocity_rad_s = aa_W_A_B.axis() * angle * inverse_dt_sec;
  Eigen::Vector3d velocity_m_s = (T_W_B.getPosition() - T_W_A.getPosition()) * inverse_dt_sec;
  // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
  DerivativeType rVal;
  rVal << velocity_m_s, angularVelocity_rad_s;
  return rVal;
}
void CubicHermiteSE3Curve::extend(const std::vector<Time>& times,
                                  const std::vector<ValueType>& values,
                                  std::vector<Key>* outKeys) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";

  // construct the Hemrite coefficients
  std::vector<Coefficient> coefficients;
  DerivativeType derivative;
  //cases:
  for (size_t i = 0; i < times.size(); ++i) {

    // only 1 new value in extend
    if (times.size() == 1) {
      // manager is empty
      if (manager_.size() == 0) {
        derivative << 0,0,0,0,0,0;

        // 1 value in manager (2 in total)
      } else if (manager_.size() == 1) {
        // get latest coefficient from manager
        CoefficientIter last = manager_.coefficientBegin();
        Time t = last->first;
        Time t1 = times[0];
        SE3 se = last->second.coefficient.getTransformation();
        SE3 se1 = values[0];
        // calculate slope
        // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
        derivative = calculateSlope(last->first,
                                    times[0],
                                    last->second.coefficient.getTransformation(),
                                    values[0]);
        // update previous coefficient
        Coefficient updated(last->second.coefficient.getTransformation(), derivative);
        manager_.updateCoefficientByKey(last->second.key, updated);
        // more than 1 values in manager
      } else if (manager_.size() > 1) {
        // get latest 2 coefficients from manager
        CoefficientIter rVal0, rVal1;
        manager_.getCoefficientsAt(manager_.coefficientEnd()->first,
                                   &rVal0,
                                   &rVal1);

        // update derivative of previous coefficient
        DerivativeType derivative0;
        derivative0 << calculateSlope(rVal0->first,
                                      times[0],
                                      rVal0->second.coefficient.getTransformation(),
                                      values[0]);
        Coefficient updated(rVal1->second.coefficient.getTransformation(), derivative0);
        manager_.updateCoefficientByKey(rVal1->second.key, updated);

        // calculate slope
        derivative << calculateSlope(rVal1->first,
                                     times[0],
                                     rVal1->second.coefficient.getTransformation(),
                                     values[0]);
      }
      // only 2 new values in extend
    } else if (times.size() == 2) {

      // manager is empty
      if (manager_.size() == 0) {

        derivative << calculateSlope(times[0],
                                     times[1],
                                     values[0],
                                     values[1]);

        // manager has at least one value
      } else {

        if (i == 0) {
          // get latest coefficient from manager
          CoefficientIter last = manager_.coefficientEnd();
          derivative << calculateSlope(last->first,
                                       times[1],
                                       last->second.coefficient.getTransformation(),
                                       values[1]);
        } else {
          derivative << calculateSlope(times[0],
                                       times[1],
                                       values[0],
                                       values[1]);
        }
      }

      // several new values in extend
    } else if (times.size() > 2) {

      if (i == 0) {
        // it is the first coefficient in curve
        if (manager_.size() == 0) {

          derivative << calculateSlope(times[0],
                                       times[1],
                                       values[0],
                                       values[1]);
        } else {
          CoefficientIter last = manager_.coefficientEnd();
          derivative << calculateSlope(last->first,
                                       times[1],
                                       last->second.coefficient.getTransformation(),
                                       values[1]);
        }

        // i is between coefficients
      } else if (i < times.size()-1) {
        // get latest coefficient from manager
        derivative << calculateSlope(times[i-1],
                                     times[i+1],
                                     values[i-1],
                                     values[i+1]);
        // i is last coefficient
      } else {
        derivative << calculateSlope(times[i-1],
                                     times[i],
                                     values[i-1],
                                     values[i]);

      }
    }
    coefficients.push_back(Coefficient(values[i], derivative));
  }
  manager_.insertCoefficients(times, coefficients, outKeys);

}

typename CubicHermiteSE3Curve::DerivativeType
CubicHermiteSE3Curve::evaluateDerivative(Time time,
                                         unsigned derivativeOrder) const {
  CHECK(false) << "Not implemented";
}

/// \brief forms cubic Hermite interpolation into a binary expression with 2 leafs and binds alpha into it,
///        uses break down of expression into its operations
gtsam::Expression<typename CubicHermiteSE3Curve::ValueType>
CubicHermiteSE3Curve::getValueExpression(const Time& time) const {

  using namespace gtsam;

  // CoefficientType is HermiteTransformation
  CoefficientIter rval0, rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;

  // shall return the Transformation Expression for Hermite
  EHermiteTransformation leaf1(rval0->second.key);
  EHermiteTransformation leaf2(rval1->second.key);

  double dt_sec = (rval1->first - rval0->first) * 1e-9;
  double alpha = double(time - rval0->first)/(rval1->first - rval0->first);

  // construct the necessary Expressions (and underlying Jacobians) for evaluation
  ETransformation transformation1 = kindr::minimal::transformationFromHermiteTransformation(leaf1);
  EVector3 angVel1 = kindr::minimal::angularVelocitiesFromHermiteTransformation(leaf1);
  EVector3 vel1 = kindr::minimal::velocitiesFromHermiteTransformation(leaf1);
  EQuaternion quat1 = kindr::minimal::rotationFromTransformation(transformation1);
  EVector3 trans1 = kindr::minimal::translationFromTransformation(transformation1);

  ETransformation transformation2 = kindr::minimal::transformationFromHermiteTransformation(leaf2);
  EVector3 angVel2 = kindr::minimal::angularVelocitiesFromHermiteTransformation(leaf2);
  EVector3 vel2 = kindr::minimal::velocitiesFromHermiteTransformation(leaf2);
  EQuaternion quat2 = kindr::minimal::rotationFromTransformation(transformation2);
  EVector3 trans2 = kindr::minimal::translationFromTransformation(transformation2);

  EQuaternion quat = kindr::minimal::hermiteQuaternionInterpolation(quat1,
                                                          angVel1,
                                                          quat2,
                                                          angVel2,
                                                          alpha,
                                                          dt_sec);

  EVector3 trans = kindr::minimal::hermiteInterpolation<Vector3>(trans1,
                                                                   vel1,
                                                                   trans2,
                                                                   vel2,
                                                                   alpha,
                                                                   dt_sec);

  return kindr::minimal::transformationFromComponents(quat, trans);
}

gtsam::Expression<typename CubicHermiteSE3Curve::DerivativeType>
CubicHermiteSE3Curve::getDerivativeExpression(const Time& time, unsigned derivativeOrder) const {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

SE3 CubicHermiteSE3Curve::evaluate(Time time) const {
  // Check if the curve is only defined at this one time
  if (manager_.getMaxTime() == time && manager_.getMinTime() == time) {
    return manager_.coefficientBegin()->second.coefficient.getTransformation();
  } else {
    CoefficientIter a, b;
    bool success = manager_.getCoefficientsAt(time, &a, &b);
    CHECK(success) << "Unable to get the coefficients at time " << time;

    // read out transformation from coefficient
    SE3 T_W_A = a->second.coefficient.getTransformation();
    SE3 T_W_B = b->second.coefficient.getTransformation();

    // read out derivative from coefficient
    gtsam::Vector6 d_W_A = a->second.coefficient.getTransformationDerivative();
    gtsam::Vector6 d_W_B = b->second.coefficient.getTransformationDerivative();

    // make alpha
    double dt_sec = (b->first - a->first) * 1e-9;
    double alpha = double(time - a->first)/(b->first - a->first);

    // Implemantation of Hermite Interpolation not easy and not fun (without expressions)!

    // translational part (easy):
    double alpha2 = alpha * alpha;
    double alpha3 = alpha2 * alpha;

    double beta0 = 2.0 * alpha3 - 3.0 * alpha2 + 1.0;
    double beta1 = -2.0 * alpha3 + 3.0 * alpha2;
    double beta2 = alpha3 - 2.0 * alpha2 + alpha;
    double beta3 = alpha3 - alpha2;

    Eigen::Vector3d translation = T_W_A.getPosition() * beta0 + T_W_B.getPosition() * beta1 +
        d_W_A.head<3>() * (beta2 * dt_sec) + d_W_B.head<3>() * (beta3 * dt_sec);

    // rotational part (hard):
    using namespace kindr::minimal;
    const double one_third = 1.0 / 3.0;
    Eigen::Vector3d scaled_d_W_A = vectorScalingImplementation<int(3)>(d_W_A.tail<3>(), one_third * dt_sec, boost::none, boost::none);
    Eigen::Vector3d scaled_d_W_B = vectorScalingImplementation<int(3)>(d_W_B.tail<3>(), one_third * dt_sec, boost::none, boost::none);

    Eigen::Vector3d w1 = inverse_rotate_point(T_W_A.getRotation(), scaled_d_W_A, boost::none, boost::none);
    Eigen::Vector3d w3 = inverse_rotate_point(T_W_B.getRotation(), scaled_d_W_B, boost::none, boost::none);
    RotationQuaternion inverse = invert_rotation_quaternion(T_W_A.getRotation(), boost::none);
    RotationQuaternion expW1 = rotationExpImplementation(-w1, boost::none);
    RotationQuaternion expW3 = rotationExpImplementation(-w3, boost::none);

    RotationQuaternion expW1_Inv = compose_rotation_quaternion(expW1, inverse, boost::none, boost::none);
    RotationQuaternion expW1_Inv_qWB = compose_rotation_quaternion(expW1_Inv, T_W_B.getRotation(), boost::none, boost::none);
    RotationQuaternion expW1_Inv_qWB_expW3 = compose_rotation_quaternion(expW1_Inv_qWB, expW3, boost::none, boost::none);

    Eigen::Vector3d w2 = rotationLogImplementation(expW1_Inv_qWB_expW3, boost::none);

    double dBeta1 = alpha3 - 3.0 * alpha2 + 3.0 * alpha;
    double dBeta2 = -2.0 * alpha3 + 3.0 * alpha2;
    double dBeta3 = alpha3;

    SO3 w1_dBeta1_exp = rotationExpImplementation(vectorScalingImplementation<int(3)>(w1, dBeta1, boost::none, boost::none), boost::none);
    SO3 w2_dBeta2_exp = rotationExpImplementation(vectorScalingImplementation<int(3)>(w2, dBeta2, boost::none, boost::none), boost::none);
    SO3 w3_dBeta3_exp = rotationExpImplementation(vectorScalingImplementation<int(3)>(w3, dBeta3, boost::none, boost::none), boost::none);

    RotationQuaternion rotation = compose_rotation_quaternion(
        compose_rotation_quaternion(
            compose_rotation_quaternion(
                T_W_A.getRotation(), w1_dBeta1_exp, boost::none, boost::none),
                w2_dBeta2_exp, boost::none, boost::none),
                w3_dBeta3_exp,boost::none, boost::none);

    return QuatTransformation(rotation, translation);
  }
}

void CubicHermiteSE3Curve::addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const {

  gtsam::noiseModel::Constrained::shared_ptr priorNoise = gtsam::noiseModel::Constrained::All(gtsam::traits<Coefficient>::dimension);

  //Add one fixed prior at priorTime and two before to ensure that at least two
  CoefficientIter rVal0, rVal1;
  manager_.getCoefficientsAt(priorTime, &rVal0, &rVal1);

  gtsam::ExpressionFactor<Coefficient> f0(priorNoise,
                                          rVal0->second.coefficient,
                                          gtsam::Expression<Coefficient>(rVal0->second.key));
  gtsam::ExpressionFactor<Coefficient> f1(priorNoise,
                                          rVal1->second.coefficient,
                                          gtsam::Expression<Coefficient>(rVal1->second.key));
  graph->push_back(f0);
  graph->push_back(f1);

}

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

void CubicHermiteSE3Curve::initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) const {
  manager_.initializeGTSAMValues(keys, values);
}

void CubicHermiteSE3Curve::initializeGTSAMValues(gtsam::Values* values) const {
  manager_.initializeGTSAMValues(values);
}

void CubicHermiteSE3Curve::updateFromGTSAMValues(const gtsam::Values& values) {
  manager_.updateFromGTSAMValues(values);
}

void CubicHermiteSE3Curve::clear() {
  manager_.clear();
}

} // namespace curves
