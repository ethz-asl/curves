/*
 * CubicHermiteSE3CurveTest.cpp
 *
 *  Created on: May 27, 2015
 *      Author: Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <gtest/gtest.h>

#include "curves/CubicHermiteSE3Curve.hpp"

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

TEST(InvarianceUnderCoordiateTransformation, Translation)
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

TEST(Debugging, FreeGaitTorsoControl)
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
