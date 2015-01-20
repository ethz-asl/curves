/*
 * @file test_Expressions_Slerp_SE3.cpp
 * @date Nov 05, 2014
 * @author Abel Gawel, Renaud Dube
 */

#include <gtest/gtest.h>
#include <curves/SlerpSE3Curve.hpp>
#include "gtsam/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "gtsam/nonlinear/ExpressionFactor.h"
#include <gtsam/base/numericalDerivative.h>

#include <Eigen/Core>
#include <boost/assign/list_of.hpp>

#define DIM 6

using namespace curves;
using namespace gtsam;
using namespace std;

typedef typename SlerpSE3Curve::ValueType ValueType;
typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;


// Expression to calculate relative measurement between 2 SE3 values
// \todo AG move this to trajectory maintainer at some point (this replaces relative pose factor functionality)
ValueType relativeMeasurementExpression(const ValueType& interp1,
                                        const ValueType& interp2,
                                        gtsam::OptionalJacobian<6,6> H1,
                                        gtsam::OptionalJacobian<6,6> H2) {
  if (H1) { *H1 = - Eigen::Matrix<double,6,6>::Identity(); }
  if (H2) { *H2 = Eigen::Matrix<double,6,6>::Identity(); }

  return ValueType(gtsam::traits<ValueType>::Local(interp1, interp2));
}

// test for correct keys and evaluation function in Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3Evaluation) {
  SlerpSE3Curve curve;
  const double t[] = {0, 10};
  const double evalTime = 5;

  // setup two SE3 objects
  ValueType poseA(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  ValueType poseB(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(2,2,2));

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
  // read out SE3 object from values container
  ValueType result = expression.value(gtsamValues);
  Eigen::Vector3d resultPos = result.getPosition();
  Eigen::Vector4d resultRot = result.getRotation().vector();

  // assert return values are as expected
  ASSERT_EQ(resultPos, Eigen::Vector3d(1,1,1));
  ASSERT_NEAR(resultRot(0),0.9238795325112867,1e-6);
  ASSERT_NEAR(resultRot(1),0.0,1e-6);
  ASSERT_NEAR(resultRot(2),-0.3826834323650897,1e-6);
  ASSERT_NEAR(resultRot(3),0.0,1e-6);
}

// test for correct keys and evaluation function in Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3CurveEvaluate) {
  SlerpSE3Curve curve;
  const double t[] = {0, 10};

  // setup two SE3 objects
  ValueType poseA(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  ValueType poseB(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(2,2,2));

  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(poseA);
  values.push_back(poseB);

  // interpolate curve
  curve.fitCurve(times, values);

  // read out SE3 object from values container
  ValueType result = curve.evaluate(5);
  Eigen::Vector3d resultPos = result.getPosition();
  Eigen::Vector4d resultRot = result.getRotation().vector();

  // assert return values are as expected
  ASSERT_EQ(resultPos, Eigen::Vector3d(1,1,1));
  ASSERT_NEAR(resultRot(0),0.9238795325112867,1e-6);
  ASSERT_NEAR(resultRot(1),0.0,1e-6);
  ASSERT_NEAR(resultRot(2),-0.3826834323650897,1e-6);
  ASSERT_NEAR(resultRot(3),0.0,1e-6);

  result = curve.evaluate(0);
  resultPos = result.getPosition();
  resultRot = result.getRotation().vector();

  ASSERT_EQ(resultPos, Eigen::Vector3d(0,0,0));
  ASSERT_NEAR(resultRot(0),1,1e-6);
  ASSERT_NEAR(resultRot(1),0,1e-6);
  ASSERT_NEAR(resultRot(2),0,1e-6);
  ASSERT_NEAR(resultRot(3),0,1e-6);

  result = curve.evaluate(10);
  resultPos = result.getPosition();
  resultRot = result.getRotation().vector();

  ASSERT_EQ(resultPos, Eigen::Vector3d(2,2,2));
  ASSERT_NEAR(resultRot(0),0.7071067811865476,1e-6);
  ASSERT_NEAR(resultRot(1),0,1e-6);
  ASSERT_NEAR(resultRot(2),-0.7071067811865476,1e-6);
  ASSERT_NEAR(resultRot(3),0,1e-6);

}

// test basic gtsam interface of Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3ExpressionGTSAMoptimization) {

  SlerpSE3Curve curve;

  // Populate coefficients
  const double t[] = {0, 45, 90};
  std::vector<Time> times(t,t+3);
  std::vector<ValueType> coefficients;
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  coefficients.push_back(ValueType(SO3(0.9238795325112867,SO3::Vector3(0,-0.3826834323650897,0)),SE3::Position(4.5,4.5,4.5)));
  coefficients.push_back(ValueType(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(9,9,9)));

  // Populate measurements
  const double tmeas[] = {20, 40, 60, 80};
  std::vector<Time> measTimes(tmeas,tmeas+4);
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(SO3(0.984807753012208, SO3::Vector3(0,-0.17364817766693036,0)),SE3::Position(2,2,2)));
  measurements.push_back(ValueType(SO3(0.9396926207859083,SO3::Vector3(0,-0.3420201433256687,0)),SE3::Position(4,4,4)));
  measurements.push_back(ValueType(SO3(0.8660254037844386,SO3::Vector3(0,-0.5,0)),SE3::Position(6,6,6)));
  measurements.push_back(ValueType(SO3(0.7660444431189781,SO3::Vector3(0,-0.6427876096865393,0)),SE3::Position(8,8,8)));

  // Fit a curve
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(times, coefficients, &outKeys);

  // Populate GTSAM values
  Values initials, expected;
  curve.initializeGTSAMValues(&initials);
  for (size_t i = 0; i < outKeys.size(); ++i) {
    expected.insert(outKeys[i],coefficients[i]);
  }

  // Noise models
  gtsam::Vector6 measNoise;
  measNoise << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
  gtsam::Vector6 priorNoise;
  priorNoise << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(measNoise);
  gtsam::noiseModel::Diagonal::shared_ptr priorNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(priorNoise);

  // Get expressions and build the graph
  NonlinearFactorGraph graph;
  for(size_t i=0; i < measurements.size(); i++) {

    Expression<ValueType> predicted(curve.getValueExpression(measTimes[i]));

    ExpressionFactor<ValueType> f(measNoiseModel, measurements[i], predicted);
    graph.add(f);

    // Assert that error is null for expected values
    std::vector<gtsam::Matrix> H(2);
    Vector error = f.unwhitenedError(expected);
    CHECK((error).cwiseAbs().maxCoeff() < 1e-6);
  }
  // Add a prior of 0 on the first coefficient
  Expression<ValueType> prior(outKeys[0]);
  graph.add(ExpressionFactor<ValueType>(priorNoiseModel,ValueType(SO3(1.0,SO3::Vector3(0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0)),prior));

  // Optimize
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials).optimize();
  ASSERT_TRUE(expected.equals(result,1e-6));
}

// test for dynamic handling of duplicate keys in gtsam::ExpressionFactor
TEST(CurvesTestSuite, testSlerpSE3ExpressionDynamicKeys){
  SlerpSE3Curve curve;

  // Populate coefficients
  const double t[] = {0, 45, 90, 135};
  std::vector<Time> times(t,t+4);
  std::vector<ValueType> coefficients;
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  coefficients.push_back(ValueType(SO3(0.9238795325112867,SO3::Vector3(0,-0.3826834323650897,0)),SE3::Position(4.5,4.5,4.5)));
  coefficients.push_back(ValueType(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(9,9,9)));
  coefficients.push_back(ValueType(SO3(0.38268343236508984,SO3::Vector3(0,-0.9238795325112866,0)),SE3::Position(13.5,13.5,13.5)));

  // Populate measurements
  const double tmeas[] = {10, 30, 20};
  std::vector<Time> measTimes(tmeas,tmeas+3);
  const double durations[] = {30, 30, 100};
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(SO3(0.9659258262890683,SO3::Vector3(0,-0.25881904510252074,0)),SE3::Position(3.633974596215561,3,2.633974596215561)));
  measurements.push_back(ValueType(SO3(0.9659258262890683,SO3::Vector3(0,-0.25881904510252074,0)),SE3::Position(4.901923788646684,3,1.901923788646684)));
  measurements.push_back(ValueType(SO3(0.6427876096865394,SO3::Vector3(0,-0.766044443118978,0)),SE3::Position(14.316911861358276,10,10.377680849309444)));

  // Fit a curve
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(times, coefficients, &outKeys);

  gtsam::Vector6 measNoise;
  measNoise << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::Sigmas(measNoise);

  // Populate GTSAM values
  Values expected;
  curve.initializeGTSAMValues(&expected);

  // create ExpressionFactor which represents relative pose relation (former Relative pose factor)
  for (int i=0; i<3; ++i) {
    Expression<ValueType> relative(relativeMeasurementExpression, curve.getValueExpression(tmeas[i]), curve.getValueExpression(tmeas[i]+durations[i]));
    ExpressionFactor<ValueType> factor(measNoiseModel,measurements[i], relative);

    Vector error = factor.unwhitenedError(expected);
    CHECK((error).cwiseAbs().maxCoeff() < 1e-6);
  }
}

ExpressionFactor<ValueType>
getFactorRelativeMeasurement(const SlerpSE3Curve& curve,
                             Time timeA, Time timeB,
                             ValueType measurement,
                             noiseModel::Diagonal::shared_ptr noiseModel) {
  Expression<ValueType> TA(curve.getValueExpression(timeA));
  Expression<ValueType> TB(curve.getValueExpression(timeB));
  Expression<ValueType> predicted = kindr::minimal::invertAndCompose(TA,TB);
  return ExpressionFactor<ValueType>(noiseModel, measurement, predicted);
}

TEST(CurvesTestSuite, testSlerpSE3RelativeExpression){
  SlerpSE3Curve curve;
  // Data taken from the MITb dataset
  // Populate expected
  const double t[] = {0, 1, 2};
  std::vector<Time> times(t,t+3);
  std::vector<ValueType> expected;
  expected.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  expected.push_back(ValueType(SO3(0.999967,SO3::Vector3(0,0,0.00760239)),SE3::Position(2.03998,0.00824351,0)));
  expected.push_back(ValueType(SO3(0.999803,SO3::Vector3(0,0,0.019962)),SE3::Position(4.26771,0.0576598,0)));

  // Populate initials
  std::vector<ValueType> initials;
  initials.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  initials.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  initials.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));

  // Populate measurements
  const double tmeas[] = {1, 2};
  std::vector<Time> measTimes(tmeas,tmeas+2);
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(SO3(0.99997,SO3::Vector3(0,0,0.0072259)),SE3::Position(2.0393,0.003006,0)));
  measurements.push_back(ValueType(SO3(0.99991,SO3::Vector3(0,0,0.013399)),SE3::Position(2.2275,0.023206,0)));

  // Fit curve
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(times, initials, &outKeys);

  //Noise models
  Vector6 measNoise;
  measNoise << 0.1, 0.1, 0.001, 0.1*M_PI/180, 0.1*M_PI/180, 0.5*M_PI/180;
  noiseModel::Diagonal::shared_ptr measNoiseModel = noiseModel::Diagonal::Sigmas(measNoise);
  Vector6 priNoise;
  priNoise << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  noiseModel::Diagonal::shared_ptr priorNoise = noiseModel::Diagonal::Sigmas(priNoise);

  NonlinearFactorGraph graph;
  Values gtsamInitial, gtsamExpected;
  curve.initializeGTSAMValues(&gtsamInitial);
  for(size_t i=0; i< outKeys.size(); i++) {
    gtsamExpected.insert(outKeys[i],expected[i]);
  }

  //prior
  Expression<ValueType> predictedPrior = curve.getValueExpression(times[0]);
  ExpressionFactor<ValueType> f(priorNoise, initials[0], predictedPrior);
  graph.push_back(f);

  //relative measurements
  graph.push_back(getFactorRelativeMeasurement(curve, measTimes[0]-1, measTimes[0],
                                               measurements[0], measNoiseModel));
  graph.push_back(getFactorRelativeMeasurement(curve, measTimes[1]-1, measTimes[1],
                                               measurements[1], measNoiseModel));

  // optimize the trajectory
  gtsam::LevenbergMarquardtParams params;
  //params.setVerbosity("ERROR");
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, gtsamInitial, params).optimize();


  ASSERT_TRUE(gtsamExpected.equals(result, 0.5));
}

TEST(CurvesTestSuite, testSlerpSE3ExpressionJacobians){
  // todo AG-RD
  ASSERT_TRUE(true);
}
