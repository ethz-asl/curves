/*
 * CubicHermiteSE3Curve.cpp
 *
 *  Created on: Feb 10, 2015
 *      Author: Abel Gawel, Renaud Dube, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <curves/CubicHermiteSE3Curve.hpp>
#include <curves/SlerpSE3Curve.hpp>
#include <iostream>

namespace curves {

CubicHermiteSE3Curve::CubicHermiteSE3Curve() : SE3Curve() {
//  hermitePolicy_.setMinimumMeasurements(4);
}

CubicHermiteSE3Curve::~CubicHermiteSE3Curve() {}

void CubicHermiteSE3Curve::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "======= Cubic Hermite SE3 CURVE =========" << std::endl;
  std::cout << str << std::endl;
  std::cout << "num of coefficients: " << manager_.size() << std::endl;
  std::cout << "dimension: " << 6 << std::endl;
  std::vector<Key> keys;
  std::vector<Time> times;
  manager_.getTimes(&times);
  manager_.getKeys(&keys);
  std::cout << "curve defined between times: " << manager_.getMinTime() <<
      " and " << manager_.getMaxTime() <<std::endl;
  std::cout <<"=========================================" <<std::endl;
  for (size_t i = 0; i < manager_.size(); i++) {
    std::cout << "coefficient " << keys[i] << ": ";
    std::cout << manager_.getCoefficientByKey(keys[i]).getTransformation() << std::endl;
//    gtsam::traits<Coefficient>::Print(manager_.getCoefficientByKey(keys[i]),ss.str());
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
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

//typedef gtsam::Expression<kindr::minimal::RotationQuaternion> EQuaternion;
//typedef gtsam::Expression<Eigen::Vector3d> EVector3;
//typedef gtsam::Expression<kindr::minimal::QuatTransformation> ETransformation;
//typedef gtsam::Expression<kindr::minimal::HermiteTransformation<double>> EHermiteTransformation;

void CubicHermiteSE3Curve::fitCurve(const std::vector<Time>& times,
                                    const std::vector<ValueType>& values, std::vector<Key>* outKeys)
{
  if(times.size() != values.size())  std::cerr << "Not the same for times and values."  << std::endl;

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
        derivative << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
      } else {
        // set velocities == 0 for start point if only one coefficient
        derivative << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
      }
    } else if (i == times.size() - 1) {
      // Last key.
//      derivative = calculateSlope(times[i-1], times[i], values[i-1], values[i]);
      derivative << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    } else {
      // Other keys.
      derivative = calculateSlope(times[i-1], times[i+1], values[i-1], values[i+1]);
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
  AngleAxis aa_W_A_B(R_W_A_B);
  double angle = aa_W_A_B.angle();
  Eigen::Vector3d angularVelocity_rad_s = aa_W_A_B.axis() * angle * inverse_dt_sec;
  Eigen::Vector3d velocity_m_s = (T_W_B.getPosition().vector() - T_W_A.getPosition().vector()) * inverse_dt_sec;
  // note: unit of derivative is m/s for first 3 and rad/s for last 3 entries
  DerivativeType rVal;
  rVal << velocity_m_s, angularVelocity_rad_s;
  return rVal;
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
}

typename CubicHermiteSE3Curve::DerivativeType CubicHermiteSE3Curve::evaluateDerivative(
    Time time, unsigned derivativeOrder) const
{
  std::cerr << "Not implemented" << std::endl;
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

SE3 CubicHermiteSE3Curve::evaluate(Time time) const {
  // Check if the curve is only defined at this one time
  if (manager_.getMaxTime() == time && manager_.getMinTime() == time) {
    return manager_.coefficientBegin()->second.coefficient.getTransformation();
  } else {
    CoefficientIter a, b;
    bool success = manager_.getCoefficientsAt(time, &a, &b);
    if(!success) std::cerr << "Unable to get the coefficients at time " << time << std::endl;

    // read out transformation from coefficient
    SE3 T_W_A = a->second.coefficient.getTransformation();
    SE3 T_W_B = b->second.coefficient.getTransformation();

    // read out derivative from coefficient
    Vector6 d_W_A = a->second.coefficient.getTransformationDerivative();
    Vector6 d_W_B = b->second.coefficient.getTransformationDerivative();

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

    SE3::Position translation(T_W_A.getPosition().vector() * beta0 + T_W_B.getPosition().vector() * beta1 +
        d_W_A.head<3>() * (beta2 * dt_sec) + d_W_B.head<3>() * (beta3 * dt_sec));

    // rotational part (hard):
    const double one_third = 1.0 / 3.0;
    Eigen::Vector3d scaled_d_W_A = one_third * dt_sec * d_W_A.tail<3>();
    Eigen::Vector3d scaled_d_W_B = one_third * dt_sec * d_W_B.tail<3>();

    Eigen::Vector3d w1 = T_W_A.getRotation().inverseRotate(scaled_d_W_A);
    Eigen::Vector3d w3 = T_W_B.getRotation().inverseRotate(scaled_d_W_B);
    RotationQuaternion inverse = T_W_A.getRotation().inverted();
    RotationQuaternion expW1 = RotationQuaternion().exponentialMap(-w1);
    RotationQuaternion expW3 = RotationQuaternion().exponentialMap(-w3);

    RotationQuaternion expW1_Inv = expW1 * inverse;
    RotationQuaternion expW1_Inv_qWB = expW1_Inv * T_W_B.getRotation();
    RotationQuaternion expW1_Inv_qWB_expW3 = expW1_Inv_qWB * expW3;

    Eigen::Vector3d w2 = expW1_Inv_qWB_expW3.logarithmicMap();

    double dBeta1 = alpha3 - 3.0 * alpha2 + 3.0 * alpha;
    double dBeta2 = -2.0 * alpha3 + 3.0 * alpha2;
    double dBeta3 = alpha3;

    SO3 w1_dBeta1_exp = RotationQuaternion().exponentialMap(dBeta1 * w1);
    SO3 w2_dBeta2_exp = RotationQuaternion().exponentialMap(dBeta2 * w2);
    SO3 w3_dBeta3_exp = RotationQuaternion().exponentialMap(dBeta3 * w3);

    RotationQuaternion rotation = T_W_A.getRotation() * w1_dBeta1_exp * w2_dBeta2_exp * w3_dBeta3_exp;

    SE3 result(translation, rotation);
    return result;
  }
}

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
  std::cerr << "Not implemented." << std::endl;
}

/// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularVelocityA(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularVelocityB(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearVelocityA(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearVelocityB(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateTwistA(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateTwistB(Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d CubicHermiteSE3Curve::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateDerivativeA(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d CubicHermiteSE3Curve::evaluateDerivativeB(unsigned derivativeOrder, Time time) {
  std::cerr << "Not implemented." << std::endl;
}

//void CubicHermiteSE3Curve::initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) const {
//  manager_.initializeGTSAMValues(keys, values);
//}
//
//void CubicHermiteSE3Curve::initializeGTSAMValues(gtsam::Values* values) const {
//  manager_.initializeGTSAMValues(values);
//}
//
//void CubicHermiteSE3Curve::updateFromGTSAMValues(const gtsam::Values& values) {
//  manager_.updateFromGTSAMValues(values);
//}

//void CubicHermiteSE3Curve::setMinSamplingPeriod(Time time) {
//  hermitePolicy_.setMinSamplingPeriod(time);
//}

void CubicHermiteSE3Curve::clear() {
  manager_.clear();
}

} // namespace
