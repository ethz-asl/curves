/*
 * PolynomialSplineVectorSpaceCurve.cpp
 *
 *  Created on: Jun 6, 2016
 *      Author: Christian Gehring
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <gtest/gtest.h>

#include "curves/PolynomialSplineScalarCurve.hpp"

using namespace curves;

typedef typename curves::PolynomialSplineQuinticScalarCurve::ValueType ValueType; // kindr::HomogeneousTransformationPosition3RotationQuaternionD ValueType
typedef typename curves::Time Time;


curves::PolynomialSplineQuinticScalarCurve::ValueType finiteDifference(curves::PolynomialSplineQuinticScalarCurve& curve, curves::Time time) {
  curves::Time dt = 1.0e-6;
  auto f_P = curve.evaluate(time+dt);
  auto f_M = curve.evaluate(time-dt);
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

  EXPECT_NEAR(3.0, curve.evaluate(1.0), 1.0e-3) << "minValue";
  EXPECT_NEAR(5.0, curve.evaluate(4.0), 1.0e-3) << "maxValue";
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

  EXPECT_NEAR(initialValue, curve.evaluate(initialTime), 1.0e-3) << "initialValue";
  EXPECT_NEAR(finalValue, curve.evaluate(finalTime), 1.0e-3) << "finalValue";
  EXPECT_NEAR(initialFirstDerivativeValue, curve.evaluateDerivative(initialTime, 1), 1.0e-3) << "initialFirstDerivativeValue";
  EXPECT_NEAR(finalFirstDerivativeValue, curve.evaluateDerivative(finalTime, 1), 1.0e-3) << "finalFirstDerivativeValue";
  EXPECT_NEAR(initialSecondDerivativeValue, curve.evaluateDerivative(initialTime, 2), 1.0e-3) << "secondFirstDerivativeValue";
  EXPECT_NEAR(finalSecondDerivativeValue, curve.evaluateDerivative(finalTime, 2), 1.0e-3) << "secondFirstDerivativeValue";
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

  EXPECT_NEAR(midValue, curve.evaluate(midTime), 1.0e-3);
  EXPECT_NEAR(finiteDifference(curve, midTime), curve.evaluateDerivative(midTime, 1), 1.0e-3) << "maximum diff";
  EXPECT_NEAR(0.0, curve.evaluateDerivative(midTime, 1), 1.0e-3) << "maximum";
  EXPECT_NEAR(finiteDifference(curve, 1.4), curve.evaluateDerivative(1.4, 1), 1.0e-3) << "inbetween";

}
