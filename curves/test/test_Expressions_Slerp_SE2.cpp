/*
 * @file test_Expressions_Slerp_SE2.cpp
 * @date Nov 25, 2015
 * @author Renaud Dubé, Abel Gawel
 */

#include <gtest/gtest.h>
#include <curves/SlerpSE2Curve.hpp>
#include <curves/Pose2_Expressions.hpp>
#include "gtsam/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "gtsam/nonlinear/ExpressionFactor.h"
#include <gtsam/base/numericalDerivative.h>

#include <Eigen/Core>
#include <boost/assign/list_of.hpp>

#include <math.h>

#include <kindr/minimal/testing-gtsam.h>

#include <curves/KeyGenerator.hpp>

#define DIM 3

using namespace curves;
using namespace gtsam;
using namespace std;

typedef typename SlerpSE2Curve::ValueType ValueType;
typedef Pose2 SE2;
typedef Rot2 SO2;

const double tolerance = 1e-3;
const double fd_step = 1e-9;

// Check for consistency between curve.evaluate() and curve.getValueExpression()
TEST(CurvesTestSuite, testSlerpSE2EvaluationVSExpression) {
  SlerpSE2Curve curve;
  const double t[] = {0, 10};
  const double evalTime = 5;

  // setup two SE2 objects
  ValueType poseA(0,0,0);
  ValueType poseB(2,2,M_PI/8);

  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(poseA);
  values.push_back(poseB);

  // interpolate curve
  curve.fitCurve(times, values);

  // get expression at evaluation time
  Expression<ValueType> expression = curve.getValueExpression(evalTime);

  // fill retrieved coefficients in gtsam values container
  Values gtsamValues;
  curve.initializeGTSAMValues(&gtsamValues);

  // read out SE2 object from values container
  ValueType result = expression.value(gtsamValues);

  SE2 directResult = curve.evaluate(evalTime);

  // assert return values are as expected
  ASSERT_EQ(directResult.x(), result.x());
  ASSERT_EQ(directResult.y(), result.y());
  ASSERT_EQ(directResult.theta(), result.theta());
}
