/*
 * PolynomialSplineVectorSpaceCurve.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <gtest/gtest.h>

#include "curves/PolynomialSplineVectorSpaceCurve.hpp"

using namespace curves;

typedef typename curves::PolynomialSplineQuinticVector3Curve::ValueType ValueType; // kindr::HomogeneousTransformationPosition3RotationQuaternionD ValueType
typedef typename curves::Time Time;

TEST(PolynomialSplineQuinticVector3Curve, Overflow)
{
  PolynomialSplineQuinticVector3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  times.push_back(0.0);
  values.push_back(ValueType(0.0, 0.0, 0.0));
  times.push_back(4.0);
  values.push_back(ValueType(1.0, 1.0, 1.0));

  curve.fitCurve(times, values);

  ValueType value;
  ASSERT_TRUE(curve.evaluate(value, -0.1));
  EXPECT_NEAR(0.0, value[0], 1e-10);
  EXPECT_NEAR(0.0, value[1], 1e-10);
  EXPECT_NEAR(0.0, value[2], 1e-10);

  ASSERT_TRUE(curve.evaluate(value, 0.0));
  EXPECT_NEAR(0.0, value[0], 1e-10);
  EXPECT_NEAR(0.0, value[1], 1e-10);
  EXPECT_NEAR(0.0, value[2], 1e-10);

  ASSERT_TRUE(curve.evaluate(value, 4.0));
  EXPECT_NEAR(1.0, value[0], 1e-10);
  EXPECT_NEAR(1.0, value[1], 1e-10);
  EXPECT_NEAR(1.0, value[2], 1e-10);

  ASSERT_TRUE(curve.evaluate(value, 4.1));
  EXPECT_NEAR(1.0, value[0], 1e-10);
  EXPECT_NEAR(1.0, value[1], 1e-10);
  EXPECT_NEAR(1.0, value[2], 1e-10);
}


TEST(PolynomialSplineQuinticVector3Curve, Debugging)
{
  PolynomialSplineQuinticVector3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  times.push_back(0.0);
  values.push_back(ValueType(0.352492, -0.208961, 0.015573));
  times.push_back(0.5);
  values.push_back(ValueType(0.352492, -0.208961, 0.115573));
  times.push_back(1.0);
  values.push_back(ValueType(0.419831, -0.2154, 0.115573));
  times.push_back(1.5);
  values.push_back(ValueType(0.419831, -0.2154, -0.000102906));

  curve.fitCurve(times, values);

//  for (double t = times[0]; t <= times[3]; t += 0.0025) {
//    std::cout << t << " " << curve.evaluate(t).x() << " " << curve.evaluate(t).y() << " " << curve.evaluate(t).z() << std::endl;
//  }
//
//  EXPECT_EQ(ValueType::Position(), curve.evaluate(0.0).getPosition());
//  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(0.0).getRotation());
//  EXPECT_EQ(ValueType::Position(), curve.evaluate(0.5).getPosition());
//  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(0.5).getRotation());
//  EXPECT_EQ(ValueType::Position(), curve.evaluate(1.0).getPosition());
//  EXPECT_EQ(ValueType::Rotation(), curve.evaluate(1.0).getRotation());
}
