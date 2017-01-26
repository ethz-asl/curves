/*
 * PolynomialSplineVectorSpaceCurve.cpp
 *
 *  Created on: Jun 6, 2016
 *      Author: Christian Gehring, Peter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <gtest/gtest.h>

#include "curves/PolynomialSplineScalarCurve.hpp"

using namespace curves;

typedef typename curves::PolynomialSplineQuinticScalarCurve::ValueType ValueType;
typedef typename curves::PolynomialSplineQuinticScalarCurve::DerivativeType DerivativeType;
typedef typename curves::Time Time;

curves::PolynomialSplineQuinticScalarCurve::ValueType finiteDifference(curves::PolynomialSplineQuinticScalarCurve& curve, curves::Time time) {
  curves::Time dt = 1.0e-6;
  curves::PolynomialSplineQuinticScalarCurve::ValueType f_P;
  curve.evaluate(f_P, time+dt);
  curves::PolynomialSplineQuinticScalarCurve::ValueType f_M;
  curve.evaluate(f_M, time-dt);
  return (f_P - f_M)/(2.0*dt);
}

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
  ValueType value;
  ASSERT_TRUE(curve.evaluate(value, -0.1));
  EXPECT_NEAR(0.0, value, 1e-10);
  ASSERT_TRUE(curve.evaluate(value, 0.0));
  EXPECT_NEAR(0.0, value, 1e-10);

  ASSERT_TRUE(curve.evaluate(value, 4.0));
  EXPECT_NEAR(1.0, value, 1e-10);

  ASSERT_TRUE(curve.evaluate(value, 4.1));
  EXPECT_NEAR(1.0, value, 1e-10);

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

  ValueType value;

  EXPECT_NEAR(1.0, curve.getMinTime(), 1.0e-3) << "minTime";
  EXPECT_NEAR(4.0, curve.getMaxTime(), 1.0e-3) << "maxTime";
  ASSERT_TRUE(curve.evaluate(value, 1.0));
  EXPECT_NEAR(3.0, value, 1.0e-3) << "minValue";
  ASSERT_TRUE(curve.evaluate(value, 4.0));
  EXPECT_NEAR(5.0, value, 1.0e-3) << "maxValue";
}

TEST(PolynomialSplineQuinticScalarCurveTest, initialAndFinalConstraints)
{
  PolynomialSplineQuinticScalarCurve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  double initialTime = 1.0;
  double initialValue = 3.0;
  double initialFirstDerivativeValue = 0.1;
  double initialSecondDerivativeValue = 0.3;

  double finalTime = 4.0;
  double finalValue = 5.0;
  double finalFirstDerivativeValue = 0.1;
  double finalSecondDerivativeValue = 0.3;

  times.push_back(initialTime);
  values.push_back(ValueType(initialValue));
  times.push_back(finalTime);
  values.push_back(ValueType(finalValue));

  curve.fitCurve(times,
                 values,
                 initialFirstDerivativeValue,
                 initialSecondDerivativeValue,
                 finalFirstDerivativeValue,
                 finalSecondDerivativeValue);

  ValueType value;
  ASSERT_TRUE(curve.evaluate(value, initialTime));
  EXPECT_NEAR(initialValue, value, 1.0e-3) << "initialValue";
  ASSERT_TRUE(curve.evaluate(value, finalTime));
  EXPECT_NEAR(finalValue, value, 1.0e-3) << "finalValue";
  DerivativeType derivative;
  curve.evaluateDerivative(derivative, initialTime, 1);
  EXPECT_NEAR(initialFirstDerivativeValue, derivative, 1.0e-3) << "initialFirstDerivativeValue";
  curve.evaluateDerivative(derivative, finalTime, 1);
  EXPECT_NEAR(finalFirstDerivativeValue, derivative, 1.0e-3) << "finalFirstDerivativeValue";
  curve.evaluateDerivative(derivative, initialTime, 2);
  EXPECT_NEAR(initialSecondDerivativeValue, derivative, 1.0e-3) << "secondFirstDerivativeValue";
  curve.evaluateDerivative(derivative, finalTime, 2);
  EXPECT_NEAR(finalSecondDerivativeValue, derivative, 1.0e-3) << "secondFirstDerivativeValue";
}

TEST(PolynomialSplineQuinticScalarCurveTest, firstDerivative)
{
  PolynomialSplineQuinticScalarCurve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;

  double initialTime = 4.0;
  double initialValue = 3.0;
  double initialFirstDerivativeValue = 0.0;
  double initialSecondDerivativeValue = 0.0;

  double finalTime = 10.0;
  double finalValue = 3.0;
  double finalFirstDerivativeValue = 0.0;
  double finalSecondDerivativeValue = 0.0;

  double midTime = (finalTime-initialTime)/2.0 + initialTime;
  double midValue = 10.0;

  times.push_back(initialTime);
  values.push_back(ValueType(initialValue));

  times.push_back(midTime);
  values.push_back(ValueType(midValue));

  times.push_back(finalTime);
  values.push_back(ValueType(finalValue));

  curve.fitCurve(times,
                 values,
                 initialFirstDerivativeValue,
                 initialSecondDerivativeValue,
                 finalFirstDerivativeValue,
                 finalSecondDerivativeValue);

  ValueType value;
  ASSERT_TRUE(curve.evaluate(value, midTime));
  EXPECT_NEAR(midValue, value, 1.0e-3);
  DerivativeType derivative;
  curve.evaluateDerivative(derivative, midTime, 1);
  EXPECT_NEAR(finiteDifference(curve, midTime), derivative, 1.0e-3) << "maximum diff";
  curve.evaluateDerivative(derivative, midTime, 1);
  EXPECT_NEAR(0.0, derivative, 1.0e-3) << "maximum";
  curve.evaluateDerivative(derivative, 1.4, 1);
  EXPECT_NEAR(finiteDifference(curve, 1.4), derivative, 1.0e-3) << "inbetween";
}

TEST(PolynomialSplineQuinticScalarCurveTest, invarianceUnderOffset)
{
  // First curve.
  PolynomialSplineQuinticScalarCurve curve1;
  std::vector<Time> times;
  std::vector<ValueType> values1;
  times.push_back(0.0);
  values1.push_back(ValueType(0.0));
  times.push_back(0.41545);
  values1.push_back(ValueType(0.1));
  times.push_back(0.58075);
  values1.push_back(ValueType(0.0));
  curve1.fitCurve(times, values1);

  // Second curve (offset).
  PolynomialSplineQuinticScalarCurve curve2;
  std::vector<ValueType> values2;
  ValueType offset = -0.05;
  for (const auto& value : values1) values2.push_back(value + offset);
  curve2.fitCurve(times, values2);

  // Check.
  for (double time = times[0]; time <= times[3]; time += 0.01) {
    ValueType value1;
    curve1.evaluate(value1, time);
    ValueType value2;
    curve2.evaluate(value2, time);
    EXPECT_NEAR(value1, value2 - offset, 1.0e-7);
  }
}
