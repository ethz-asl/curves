/*
 * CubicHermiteSE3CurveTest.cpp
 *
 *  Created on: May 27, 2015
 *      Author: Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <gtest/gtest.h>

#include "curves/CubicHermiteSE3Curve.hpp"
#include <kindr/Core>
#include <kindr/common/gtest_eigen.hpp>
#include <limits>

typedef std::numeric_limits< double > dbl;

using namespace curves;

typedef typename curves::CubicHermiteSE3Curve::ValueType ValueType; // kindr::HomogeneousTransformationPosition3RotationQuaternionD ValueType
typedef typename curves::Time Time;

TEST(Evaluate, IdentityPoses)
{
  CubicHermiteSE3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  times.push_back(0.0);
  values.push_back(ValueType());
  times.push_back(1.0);
  values.push_back(ValueType());
  curve.fitCurve(times, values);

  EXPECT_EQ(ValueType::Position(), curve.evaluate(0.0).getPosition());
  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(0.0).getRotation());
  EXPECT_EQ(ValueType::Position(), curve.evaluate(0.5).getPosition());
  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(0.5).getRotation());
  EXPECT_EQ(ValueType::Position(), curve.evaluate(1.0).getPosition());
  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(1.0).getRotation());
}

TEST(Evaluate, TranslationOnly)
{
  CubicHermiteSE3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  times.push_back(0.0);
  values.push_back(ValueType(ValueType::Position(-10.0, 10.0, 10.0), ValueType::Rotation()));
  times.push_back(1.0);
  values.push_back(ValueType(ValueType::Position(10.0, -10.0, -10.0), ValueType::Rotation()));
  curve.fitCurve(times, values);

  EXPECT_EQ(values[0].getPosition(), curve.evaluate(0.0).getPosition());
  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(0.0).getRotation());
  EXPECT_EQ(ValueType::Position(0.0, 0.0, 0.0), curve.evaluate(0.5).getPosition());
  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(0.5).getRotation());
  EXPECT_EQ(values[1].getPosition(), curve.evaluate(1.0).getPosition());
  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(1.0).getRotation());
}

TEST(InvarianceUnderCoordinateTransformation, Translation)
{
  CubicHermiteSE3Curve curve1;
  std::vector<Time> times;
  std::vector<ValueType> values1;
  times.push_back(0.0);
  values1.push_back(ValueType(ValueType::Position(-1.0, 0.1, 10.0), ValueType::Rotation()));
  times.push_back(5.0);
  values1.push_back(ValueType(ValueType::Position(1.0, 0.2, 5.0), ValueType::Rotation()));
  curve1.fitCurve(times, values1);

  CubicHermiteSE3Curve curve2;
  ValueType::Position offset(10.0, 10.0, 10.0);
  std::vector<ValueType> values2;
  for (const auto& value : values1) {
    values2.push_back(value);
    values2.back().getPosition() += offset;
  }
  curve2.fitCurve(times, values2);

  for (double time = times[0]; time <= times[1]; time += 0.1) {
    const ValueType::Position position1 = curve1.evaluate(time).getPosition();
    const ValueType::Rotation rotation1 = curve1.evaluate(time).getRotation();
    const ValueType::Position position2 = curve2.evaluate(time).getPosition();
    const ValueType::Rotation rotation2 = curve2.evaluate(time).getRotation();
    EXPECT_NEAR(position1.x(), position2.x() - offset.x(), 1e-6);
    EXPECT_NEAR(position1.y(), position2.y() - offset.y(), 1e-6);
    EXPECT_NEAR(position1.z(), position2.z() - offset.z(), 1e-6);
    EXPECT_NEAR(rotation1.x(), rotation2.x(), 1e-6);
    EXPECT_NEAR(rotation1.y(), rotation2.y(), 1e-6);
    EXPECT_NEAR(rotation1.z(), rotation2.z(), 1e-6);
    EXPECT_NEAR(rotation1.w(), rotation2.w(), 1e-6);
  }
}

TEST(InvarianceUnderCoordinateTransformation, Rotation)
{
  CubicHermiteSE3Curve curve1;
  std::vector<Time> times;
  std::vector<ValueType> values1;
  times.push_back(0.0);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(kindr::EulerAnglesYprPD(0.0, 0.0, 0.0))));
  times.push_back(2.0);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(kindr::EulerAnglesYprPD(20.0 / 180.0 * M_PI, 0.0, 0.0))));
  times.push_back(4.0);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(kindr::EulerAnglesYprPD(-20.0 / 180.0 * M_PI, 0.0, 0.0))));
  times.push_back(5.2);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(kindr::EulerAnglesYprPD(0.0, 0.0, 0.0))));
  curve1.fitCurve(times, values1);

  CubicHermiteSE3Curve curve2;
  ValueType::Rotation transform(kindr::EulerAnglesYprPD(3.0, 0.0, 0.0));
  std::vector<ValueType> values2;
  for (const auto& value : values1) {
    values2.push_back(value);
    values2.back().getRotation() = transform * values2.back().getRotation();
//    values2.back().getRotation().setUnique(); // This results in a flip in w.
  }
  curve2.fitCurve(times, values2);

  for (double time = times[0]; time <= times[3]; time += 0.1) {
    const ValueType::Position position1 = curve1.evaluate(time).getPosition();
    const ValueType::Rotation rotation1 = curve1.evaluate(time).getRotation();
    const ValueType::Position position2 = curve2.evaluate(time).getPosition();
    const ValueType::Rotation rotation2 = transform.inverted() * curve2.evaluate(time).getRotation();

//    std::cout << 180.0 / M_PI * kindr::EulerAnglesYprPD(curve2.evaluate(time).getRotation()).getUnique().yaw() << ", ";
//    std::cout << "[" << curve2.evaluate(time).getRotation() << "]" << std::endl;
    EXPECT_NEAR(position1.x(), position2.x(), 1e-6);
    EXPECT_NEAR(position1.y(), position2.y(), 1e-6);
    EXPECT_NEAR(position1.z(), position2.z(), 1e-6);
    EXPECT_LT(rotation1.getDisparityAngle(rotation2), 1e-6) << "rot1: " <<  rotation1 << "  rot2: " << rotation2 << std::endl;
  }
//  std::cout << std::endl;
}

TEST(InvarianceUnderCoordinateTransformation, Rotation2)
{
  CubicHermiteSE3Curve curve1;
  std::vector<Time> times;
  std::vector<ValueType> values1;
  times.push_back(0.0);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(0.0652549, 0.0, 0.0, 0.997869)));
  times.push_back(2.0);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(0.109015, 0.0, 0.0, -0.99404))); // w positive.
  times.push_back(4.0);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(0.237542, 0.0, 0.0, 0.971377)));
  times.push_back(5.2);
  values1.push_back(ValueType(ValueType::Position(), ValueType::Rotation(0.0652549, 0.0, 0.0, 0.997869)));
  curve1.fitCurve(times, values1);

  CubicHermiteSE3Curve curve2;
  ValueType::Rotation transform = values1.front().getRotation();
  std::vector<ValueType> values2;
  for (const auto& value : values1) {
    values2.push_back(value);
    values2.back().getRotation() = transform * values2.back().getRotation();
//    values2.back().getRotation().setUnique(); // Enforces w positive.
  }
  curve2.fitCurve(times, values2);

  for (double time = times[0]; time <= times[3]; time += 0.1) {
    const ValueType::Position position1 = curve1.evaluate(time).getPosition();
    const ValueType::Rotation rotation1 = curve1.evaluate(time).getRotation();
    const ValueType::Position position2 = curve2.evaluate(time).getPosition();
    const ValueType::Rotation rotation2 = transform.inverted() * curve2.evaluate(time).getRotation();
//    std::cout << 180.0 / M_PI * kindr::EulerAnglesYprPD(curve1.evaluate(time).getRotation()).getUnique().yaw() << ", ";
//    std::cout << "[" <<curve1.evaluate(time).getRotation() << "]" << std::endl;
    EXPECT_NEAR(position1.x(), position2.x(), 1e-6);
    EXPECT_NEAR(position1.y(), position2.y(), 1e-6);
    EXPECT_NEAR(position1.z(), position2.z(), 1e-6);
    EXPECT_LT(rotation1.getDisparityAngle(rotation2), 1e-6);
  }
//  std::cout << std::endl;
}


TEST(CubicHermiteSE3CurveTest, firstDerivative)
{
  CubicHermiteSE3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  double time0 = -1.56;
  times.push_back(time0);
  ValueType::Rotation rotation0(kindr::EulerAnglesZyxD(M_PI_2, 0.2, -0.9));
  ValueType::Position position0(1.0, 2.0, 4.0);
  ValueType transform0(position0, rotation0);
  values.push_back(transform0);

  double timeMid1 = 1.0;
  times.push_back(timeMid1);
  ValueType::Rotation rotationMid1(kindr::EulerAnglesZyxD(2.0, 3.0, -1.1));
  ValueType::Position positionMid1(2.0, 4.0, 8.0);
  values.push_back(ValueType(positionMid1, rotationMid1));

  double timeMid2 = 2.5;
  times.push_back(timeMid2);
  ValueType::Rotation rotationMid2(kindr::EulerAnglesZyxD(0.2, 0.5, 0.2));
  ValueType::Position positionMid2(2.0, 4.0, 8.0);
  values.push_back(ValueType(positionMid2, rotationMid2));


  double timeMid3 = 3.0;
  times.push_back(timeMid3);
  ValueType::Rotation rotationMid3(kindr::EulerAnglesZyxD(1.5, 0.4, -0.3));
  ValueType::Position positionMid3(2.0, 4.0, 8.0);
  values.push_back(ValueType(positionMid3, rotationMid3));


  double time1 = 4.0;
  times.push_back(time1);
  ValueType::Rotation rotation1(kindr::EulerAnglesZyxD(0.0, 0.0, 0.0));
  ValueType::Position position1(4.0, 8.0, 16.0);
  values.push_back(ValueType(position1, rotation1.getUnique()));
  curve.fitCurve(times, values);

  // Check first knot
  ValueType transform = curve.evaluate(time0);
  ValueType expTransform = transform0;
  EXPECT_NEAR(expTransform.getPosition().x(), transform.getPosition().x(), 1e-6);
  EXPECT_NEAR(expTransform.getPosition().y(), transform.getPosition().y(), 1e-6);
  EXPECT_NEAR(expTransform.getPosition().z(), transform.getPosition().z(), 1e-6);
  EXPECT_NEAR(0.0, expTransform.getRotation().getDisparityAngle(transform.getRotation()), 1e-3);

  // Derivative at first knot
  CubicHermiteSE3Curve::DerivativeType derivative = curve.evaluateDerivative(time0, 1);
  CubicHermiteSE3Curve::DerivativeType expDerivative;
  expDerivative << Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero();
  KINDR_ASSERT_DOUBLE_MX_EQ(expDerivative, derivative, 1e-1, "first");

  // Derivative at last knot
  derivative = curve.evaluateDerivative(time1, 1);
  expDerivative << Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero();
  KINDR_ASSERT_DOUBLE_MX_EQ(expDerivative, derivative, 1e-1, "last");


  // Finite difference
  double h = 1.0e-8;

  for (double time = times[0]+h; time <= times[3]; time += 0.1) {
//    double time = 2.8;//0.7;
    // double time = timeMid+1.0e-8; does not work 0.0 1e0 2e0  ;
    double timeA = time-h;
    double timeB = time+h;
    ValueType T_A = curve.evaluate(timeA);
    ValueType T_B = curve.evaluate(timeB);
    Eigen::Vector3d angularVel = T_B.getRotation().getUnique().boxMinus(T_A.getRotation().getUnique())/(2.0*h);
    Eigen::Vector3d angularVel2 = T_B.getRotation().boxMinus(T_A.getRotation())/(2.0*h);

  //  std::cout << "+===================================================" << angularVel.transpose() << std::endl;
  //  std::cout << "angularVel: " << angularVel.transpose() << std::endl;
  //  std::cout << "angularVel2: " << angularVel2.transpose() << std::endl;
  //  std::cout.precision(dbl::max_digits10);
  //  std::cout << "rotA " << T_A.getRotation() << std::endl;
  //  std::cout << "rotA u " << T_A.getRotation().getUnique() << std::endl;
  //  std::cout << "rotB " << T_B.getRotation() << std::endl;
  //  std::cout << "rotB u " << T_B.getRotation().getUnique() << std::endl;
  //  std::cout << "norm: " << angularVel.norm() << " pi./2 " << M_PI_2 << std::endl;
    Eigen::Vector3d linearVel = (T_B.getPosition().vector() - T_A.getPosition().vector())/(2.0*h);
  //  std::cout << "======> test\n";
  //  std::cout << "T_A: " << T_A << std::endl;
  //  std::cout << "T_B: " << T_B << std::endl;
  //  std::cout << "linearVel: " << linearVel << std::endl;

    expDerivative << linearVel, angularVel;
    derivative = curve.evaluateDerivative(time, 1);
    KINDR_ASSERT_DOUBLE_MX_EQ_ZT(expDerivative, derivative, 1.0, "fd", 1.0e-7);
  //  EXPECT_NEAR(expDerivative(0), derivative(0), 1e-2);
  //  EXPECT_NEAR(expDerivative(1), derivative(1), 1e-2);
  //  EXPECT_NEAR(expDerivative(2), derivative(2), 1e-2);
  //  EXPECT_NEAR(expDerivative(3), derivative(3), 1e-2);
  //  EXPECT_NEAR(expDerivative(4), derivative(4), 1e-2);
  //  EXPECT_NEAR(expDerivative(5), derivative(5), 1e-2);
  //  std::cout << "derivative: " << derivative << std::endl;
  //  std::cout << "expexpDerivative: " << expDerivative << std::endl;
  }
}

TEST(Debugging, DISABLED_FreeGaitTorsoControl)
{
  CubicHermiteSE3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;
  ValueType::Rotation fixRotation(0.999999, -6.31036e-05, 0.00109732, -3.43683e-06);
//  fixRotation.setIdentity(); // With this it succeeds.
  times.push_back(0.0);
  values.push_back(ValueType(ValueType::Position(99.9, 0.0, 0.38), fixRotation));
  times.push_back(1.5);
  values.push_back(ValueType(ValueType::Position(100.0, 0.0, 0.40), fixRotation));
  curve.fitCurve(times, values);

//  curve.print();

  for (double time = times[0]; time <= times[1]; time += 0.1) {

    const ValueType::Position position = curve.evaluate(time).getPosition();
    const ValueType::Rotation rotation = curve.evaluate(time).getRotation();

//    std::cout << "Time: " << time << " s, Position: " << position << std::endl;
//    std::cout << "Time: " << time << " s, Rotation: " << rotation << std::endl;

    EXPECT_GE(position.x(), values[0].getPosition().x());
//    EXPECT_GE(position.y(), values[0].getPosition().y()); // Numerical issues.
    EXPECT_GE(position.z(), values[0].getPosition().z());
    EXPECT_LE(position.x(), values[1].getPosition().x());
//    EXPECT_LE(position.y(), values[1].getPosition().y()); // Numerical issues.
    EXPECT_LE(position.z(), values[1].getPosition().z());

    EXPECT_NEAR(fixRotation.x(), rotation.x(), 1e-6);
    EXPECT_NEAR(fixRotation.y(), rotation.y(), 1e-6);
    EXPECT_NEAR(fixRotation.z(), rotation.z(), 1e-6);
    EXPECT_NEAR(fixRotation.w(), rotation.w(), 1e-6);
  }
}

TEST(GetTime, Simple)
{
  CubicHermiteSE3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  times.push_back(1.1);
  values.push_back(ValueType());
  times.push_back(5.6);
  values.push_back(ValueType());
  times.push_back(7.2);
  values.push_back(ValueType());
  curve.fitCurve(times, values);

  EXPECT_EQ(times[0], curve.getMinTime());
  EXPECT_EQ(times[2], curve.getMaxTime());
}
