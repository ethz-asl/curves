/*
 * @file test_Expressions.cpp
 * @date Oct 31, 2014
 * @author Renaud Dube, Abel Gawel
 */

#include <gtest/gtest.h>
#include <curves/LinearInterpolationVectorSpaceCurve.hpp>
#include <gtsam/nonlinear/Expression.h>
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/ExpressionFactor.h>

#include <boost/assign/list_of.hpp>

#define DIM 3

using namespace curves;
using namespace gtsam;
using namespace std;

typedef typename LinearInterpolationVectorSpaceCurve<DIM>::ValueType ValueType;

// test correct handling of keys and interpolation with expressions
TEST(CurvesTestSuite, testExpressionKeysAndEvaluation) {

  // create curve
  LinearInterpolationVectorSpaceCurve<DIM> curve;
  const double t[] = {0, 10};
  const double evalTime = 5;
  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(ValueType(0,0,0));
  values.push_back(ValueType(2,2,2));
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(times, values, &outKeys);
  ASSERT_EQ(outKeys.size(), values.size());

  // Retrieve Expression from curve
  Expression<ValueType> expression = curve.getEvalExpression(evalTime);

  // populate gtsam Values container with coefficients & keys, 3 ways:
  Values gtsamValues;
  gtsam::FastVector<gtsam::Key> keys;

  // classic approach
  for (int i=0; i<values.size(); ++i) {
    gtsamValues.insert(outKeys[i], ValueType(values[i]));
    keys.push_back(outKeys[i]);
  }
  gtsamValues.clear();

  // populate all values
  curve.initializeGTSAMValues(&gtsamValues);
  gtsamValues.clear();

  //populate selected values
  curve.initializeGTSAMValues(keys, &gtsamValues);

  // read out interpolated value from gtsam Values container via expression
  ValueType result = expression.value(gtsamValues);
  // Check evaluation results
  ASSERT_EQ(result, Eigen::Vector3d(1,1,1));
}

// test optimization of gtsam factor graph using expressions
TEST(CurvesTestSuite, testExpressionGTSAMoptimization) {

  LinearInterpolationVectorSpaceCurve<DIM> curve;

  // Populate coefficients
  const double t[] = {0, 10, 20};
  const double evalTime = 5;
  std::vector<Time> times(t,t+3);
  std::vector<ValueType> coefficients;
  coefficients.push_back(ValueType(0,0,0));
  coefficients.push_back(ValueType(40,20,-30));
  coefficients.push_back(ValueType(20,-20,-80));

  // Populate measurements
  const double tmeas[] = {4, 8, 12, 16};
  std::vector<Time> measTimes(tmeas,tmeas+4);
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(16,8,-12));
  measurements.push_back(ValueType(32,16,-24));
  measurements.push_back(ValueType(36,12,-40));
  measurements.push_back(ValueType(28,-4,-60));

  // Fit a curve
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(times, coefficients, &outKeys);

  // Populate GTSAM values
  Values initials, expected;
  curve.initializeGTSAMValues(outKeys, &expected);
  for(int i=0; i< coefficients.size(); i++) {
    initials.insert(outKeys[i],ValueType(0,0,0));
  }

  // Noise models
  SharedNoiseModel measNoiseModel = noiseModel::Diagonal::Sigmas(ValueType(1, 1, 1));
  SharedNoiseModel priorNoiseModel = noiseModel::Diagonal::Sigmas(ValueType(0, 0, 0));

  // Get expressions and build the graph
  NonlinearFactorGraph graph;
  for(int i=0; i < measurements.size(); i++) {
    Expression<ValueType> predicted(curve.getEvalExpression(measTimes[i]));

    ExpressionFactor<ValueType > f(measNoiseModel,
                                   ValueType(measurements[i]),
                                   predicted);
    graph.add(f);
    // Assert that error is null for expected values
    std::vector<Matrix> H(2);
    Vector error = f.unwhitenedError(expected, H);
    ASSERT_EQ(error, ValueType(0,0,0));
  }

  // Add a prior of 0,0,0 on the first coefficient
  Expression<ValueType> leaf1(outKeys[0]);
  Expression<ValueType > prior(leaf1);
  graph.add(ExpressionFactor<ValueType>(priorNoiseModel,ValueType(ValueType(0,0,0)),prior));

  // Optimize
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials).optimize();

  ASSERT_TRUE(expected.equals(result,1e-8));
}



