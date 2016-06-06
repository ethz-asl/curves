/*
 * PolynomialSplineVectorSpaceCurve.cpp
 *
 *  Created on: Jun 16, 2015
 *      Author: Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <gtest/gtest.h>

#include "curves/PolynomialSplineScalarCurve.hpp"

using namespace curves;

typedef typename curves::PolynomialSplineQuinticScalarCurve::ValueType ValueType; // kindr::poses::eigen_impl::HomogeneousTransformationPosition3RotationQuaternionD ValueType
typedef typename curves::Time Time;

TEST(PolynomialSplineQuinticScalarCurveTest, Overflow)
{
  PolynomialSplineQuinticScalarCurve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  times.push_back(0.0);
  values.push_back(ValueType(0.0));
  times.push_back(4.0);
  values.push_back(ValueType(1.0));

  curve.fitCurve(times, values);

  EXPECT_NEAR(0.0, curve.evaluate(-0.1), 1e-10);
  EXPECT_NEAR(0.0, curve.evaluate(0.0), 1e-10);
  EXPECT_NEAR(1.0, curve.evaluate(4.0), 1e-10);
  EXPECT_NEAR(1.0, curve.evaluate(4.1), 1e-10);
}

TEST(PolynomialSplineQuinticScalarCurveTest, minMax)
{
  PolynomialSplineQuinticScalarCurve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  times.push_back(1.0);
  values.push_back(ValueType(3.0));
  times.push_back(4.0);
  values.push_back(ValueType(5.0));

  curve.fitCurve(times, values);

  EXPECT_NEAR(1.0, curve.getMinTime(), 1.0e-3) << "minTime";
  EXPECT_NEAR(4.0, curve.getMaxTime(), 1.0e-3) << "maxTime";

  EXPECT_NEAR(1.0, curve.evaluate(3.0), 1.0e-3) << "minValue";
  EXPECT_NEAR(4.0, curve.evaluate(5.0), 1.0e-3) << "maxValue";
}
