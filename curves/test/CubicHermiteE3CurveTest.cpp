/*
 * CubicHermiteSE3CurveTest.cpp
 *
 *  Created on: May 27, 2015
 *      Author: Péter Fankhauser, Christian Gehring
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#include <gtest/gtest.h>

#include "curves/CubicHermiteE3Curve.hpp"
#include <kindr/Core>
#include <kindr/common/gtest_eigen.hpp>
#include <limits>

using namespace curves;

typedef typename curves::CubicHermiteE3Curve::ValueType ValueType;
typedef typename curves::CubicHermiteE3Curve::DerivativeType DerivativeType;
typedef typename curves::Time Time;

TEST(CubicHermiteE3Curve, TwoPoints)
{
  CubicHermiteE3Curve curve;
  std::vector<Time> times;
  std::vector<ValueType> values;
  std::vector<ValueType> result1;
  std::vector<ValueType> result2;
  ValueType value1;
  ValueType value2;

  //fit a curve through 2 points

  times.push_back(0.0);
  times.push_back(1.0);

  values.push_back(ValueType(0.0, 0.0, 0.0));
  values.push_back(ValueType(1.0, 1.0, 0.0));

  curve.fitCurve(times, values);

  for (double t = 0.0; t <= 1.0; t += 0.05) {
    curve.evaluate(value1, t);
    result1.push_back(value1);
    //std::cerr << t << ": " << value1.transpose() << std::endl;
  }

  //Use same instance of curve to compute spline with more points

  times.clear();
  times.push_back(0.0);
  times.push_back(0.25);
  times.push_back(0.50);
  times.push_back(0.75);
  times.push_back(1.0);

  values.clear();
  values.push_back(ValueType(0.0, 0.0, 0.0));
  values.push_back(ValueType(0.25, 0.25, 0.0));
  values.push_back(ValueType(0.50, 0.50, 0.0));
  values.push_back(ValueType(0.75, 0.75, 0.0));
  values.push_back(ValueType(1.0, 1.0, 0.0));

  curve.fitCurve(times, values);

  //Fit again original curve

  times.clear();
  times.push_back(0.0);
  times.push_back(1.0);

  values.clear();
  values.push_back(ValueType(0.0, 0.0, 0.0));
  values.push_back(ValueType(1.0, 1.0, 0.0));

  curve.fitCurve(times, values);

  for (double t = 0.0; t <= 1.0; t += 0.05) {
    curve.evaluate(value2, t);
    result2.push_back(value2);
    //std::cerr << t << ": " << value2.transpose() << std::endl;
  }

  //assert that result1 and result2 are the same (or very close)

  for (size_t i = 0; i < result1.size(); i++) {
    EXPECT_NEAR(result1[i].x(), result2[i].x(), 1e-6);
    EXPECT_NEAR(result1[i].y(), result2[i].y(), 1e-6);
    EXPECT_NEAR(result1[i].z(), result2[i].z(), 1e-6);
  }

}

