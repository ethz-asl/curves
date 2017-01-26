/*
 * @file test_Hermite.cpp
 * @date Feb 11, 2015
 * @author Abel Gawel, Renaud Dube
 */

#include <gtest/gtest.h>
#include <curves/CubicHermiteSE3Curve.hpp>
#include "gtsam/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "gtsam/nonlinear/ExpressionFactor.h"
#include <gtsam/base/numericalDerivative.h>
#include "kindr/minimal/testing-gtsam.h"

#include "curves/helpers.hpp"

#include <Eigen/Core>
#include <boost/assign/list_of.hpp>

#define DIM 6

using namespace curves;
using namespace gtsam;
using namespace std;

typedef typename CubicHermiteSE3Curve::ValueType ValueType;
typedef typename CubicHermiteSE3Curve::Coefficient Coefficient;
typedef typename Eigen::Matrix<double, 6,1> Vector6d;
typedef typename CubicHermiteSE3Curve::DerivativeType DerivativeType;
typedef typename kindr::minimal::HermiteTransformation<double> HermiteTransformation;
typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;

// Test the HermiteTransfoamtion Coefficient struct
// todo move this test together with the HermiteTransformation to minkindr_gtsam
TEST(CurvesTestSuite, testHermiteTransformationStruct) {

  using gtsam::Expression;
  gtsam::Values values;

  // Create some values
  for (int i=1; i<6; ++i) {
    SO3 Rval;
    Eigen::Vector3d tval;
    Eigen::Vector3d wval;
    Eigen::Vector3d vval;
    Vector6d dval;

    Rval.setRandom();
    tval.setRandom();
    wval.setRandom();
    vval.setRandom();
    dval << vval, wval;
    ValueType T(Rval, tval);

    HermiteTransformation HT(T, dval);
    values.insert(i, HT);
    kindr::minimal::EHermiteTransformation EHT(i);
    EXPECT_TRUE(EHT.value(values).getTransformation() == HT.getTransformation());
    EXPECT_TRUE(EHT.value(values).getTransformationDerivative() ==
        HT.getTransformationDerivative());

    const double fd_step = 1e-9;
    const double tolerance = 1e-6;
    SCOPED_TRACE("Testing Expression Jacobians.");
    gtsam::testExpressionJacobians(EHT, values, fd_step, tolerance);
  }
}

// Test the HermiteTransfoamtion Coefficient expressions
TEST(CurvesTestSuite, testHermiteTransformationExpressions) {
  using gtsam::Expression;
  gtsam::Values values;

  // Create some values
  SO3 Rval;
  Eigen::Vector3d tval;
  Eigen::Vector3d wval;
  Eigen::Vector3d vval;
  Vector6d dval;

  Rval.setRandom();
  tval.setRandom();
  wval.setRandom();
  vval.setRandom();
  dval << vval, wval;
  ValueType T(Rval, tval);

  HermiteTransformation HT(T, dval);
  values.insert(1, HT);
  values.insert(2, T);
  kindr::minimal::EHermiteTransformation EHT(1);
  kindr::minimal::ETransformation ET(2);

  kindr::minimal::ETransformation eTransform =
      transformationFromHermiteTransformation(EHT);

  kindr::minimal::EVector3 eAngVel =
      angularVelocitiesFromHermiteTransformation(EHT);

  kindr::minimal::EVector3 eVel =
      velocitiesFromHermiteTransformation(EHT);

  const double fd_step = 1e-9;
  const double tolerance = 1e-6;
  SCOPED_TRACE("Testing Expression Jacobians.");
  testExpressionJacobians(EHT, values, fd_step, tolerance);
  testExpressionJacobians(eTransform, values, fd_step, tolerance);
  testExpressionJacobians(eAngVel, values, fd_step, tolerance);
  testExpressionJacobians(eVel, values, fd_step, tolerance);

}

TEST(CurvesTestSuite, testCubicHermiteVelocityScale) {
  HermiteTransformation H1, H2;

  SE3 T1(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  SE3 T2(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(10,20,30));
  SE3 Ti(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(2.5,5,7.5));

  Vector3 v1, v2, vr1, vr2, vi, vri;
  v1 << 1,2,3;
  v2 = v1;
  vi = v1;
  vr1 << 0,0,0;
  vr2 = vr1;
  vri = vr1;

  double t1, t2, ti;
  t1 = 0;
  t2 = 10;
  ti = 2.5;

  double alpha = (ti-t1)/(t2-t1);
  Expression<SE3> eT1(1), eT2(2);
  Expression<Vector3> ev1(3), ev2(4), evr1(5), evr2(6);
  gtsam::Values values;
  values.insert(1, T1);
  values.insert(2, T2);
  values.insert(3, v1);
  values.insert(4, v2);
  values.insert(5, vr1);
  values.insert(6, vr2);


  Expression<Vector3> eInterpolated = kindr::minimal::hermiteInterpolation<Vector3>(
      kindr::minimal::translationFromTransformation(eT1),
      ev1,
      kindr::minimal::translationFromTransformation(eT2),
      ev2,
      alpha,
      t2 - t1);

  Vector3 TiExpectedTrans = eInterpolated.value(values);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Ti.getPosition(), TiExpectedTrans, 1e-8));

  Expression<SO3> erInterpolated = kindr::minimal::hermiteQuaternionInterpolation(
        kindr::minimal::rotationFromTransformation(eT1),
        evr1,
        kindr::minimal::rotationFromTransformation(eT2),
        evr2,
        alpha,
        t2 - t1);

  SO3 TiExpectedRot = erInterpolated.value(values);
  EXPECT_TRUE(EIGEN_MATRIX_NEAR(Ti.getRotation().vector(), TiExpectedRot.vector(), 1e-8));

}

// test for correct keys and evaluation function in Slerp SE3 curves
TEST(CurvesTestSuite, testCubicHermiteEvaluation) {

  // setup three SE3 objects
  const Time t[] = {0, Time(1e9), Time(2e9)};
  std::vector<Time> times(t,t+3);
  std::vector<SE3> values;
  for (size_t i = 0; i<3; ++i) {
    SO3 Rval;
    Eigen::Vector3d tval;
//    Rval = SO3(1,SO3::Vector3(0,0,0));
    Rval.setRandom();
    tval.setRandom();
    SE3 T(Rval, tval);
    values.push_back(T);
  }

  // interpolate curve
  CubicHermiteSE3Curve curve;
  curve.extend(times, values);
  curve.print("curve ");
  // fill retrieved coefficients in gtsam values container
  Values gtsamValues;
  curve.initializeGTSAMValues(&gtsamValues);

  // ensure that getValuexpression and evaluate produce same values
  for (Time ti = 0; ti < 2e9; ti+=0.1e9) {
    std::vector<std::string> pose, pose2;
    SE3 eval = curve.evaluate(ti);
    Expression<ValueType> expressionC = curve.getValueExpression(ti);
    ValueType resultC = expressionC.value(gtsamValues);

    if (!gtsam::traits<SE3>::Equals(eval, resultC, 1e-8)) {
      std::cout << "ti " << ti*1e-9 << std::endl;
      gtsam::traits<SE3>::Print(eval, "eval ");
      gtsam::traits<SE3>::Print(resultC, "expression ");
    }
    EXPECT_TRUE(gtsam::traits<SE3>::Equals(eval, resultC, 1e-8));

  }
}

#if 0
// test optimization of gtsam factor graph using expressions
TEST(CurvesTestSuite, testHermiteExpressionGTSAMoptimization) {

  CubicHermiteSE3Curve curve;

  // Populate coefficients
  const double t[] = {0, 10, 20, 30, 40};
  std::vector<Time> times(t,t+5);
  std::vector<ValueType> coefficients;
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(10,0,0)));
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(20,0,0)));
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(30,0,0)));
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(40,0,0)));

  // Populate measurements
  const double tmeas[] = {4, 8, 12, 16, 20, 24, 28, 32, 36};
  std::vector<Time> measTimes(tmeas,tmeas+9);
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(4,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(8,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(12,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(16,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(20,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(24,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(28,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(32,0,0)));
  measurements.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(36.0,0,0)));

  // Fit a curve
  std::vector<gtsam::Key> outKeys;
  curve.extend(times, coefficients, &outKeys);

  // Populate GTSAM values
  Values initials, expected;
  curve.initializeGTSAMValues(outKeys, &expected);

  for(size_t i=0; i< coefficients.size(); i++) {
    Eigen::Matrix<double, 6, 1> derivatives;
    derivatives << 0,0,0,0,0,0;
    ValueType pose(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
    kindr::minimal::HermiteTransformation<double> init(pose, derivatives);
    initials.insert(outKeys[i],init);
  }

  // Noise models
  gtsam::Vector6 measNoise;
  Eigen::Matrix<double, 6, 1> priorNoise;
  priorNoise << 0,0,0,0,0,0;
  measNoise << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
  SharedNoiseModel measNoiseModel = noiseModel::Diagonal::Sigmas(measNoise);
  SharedNoiseModel priorNoiseModel = noiseModel::Diagonal::Sigmas(priorNoise);

  // Get expressions and build the graph
  NonlinearFactorGraph graph;
  for(size_t i=0; i < measurements.size(); i++) {
    Expression<ValueType> predicted(curve.getValueExpression(measTimes[i]));

    ExpressionFactor<ValueType> f(measNoiseModel,
                                  ValueType(measurements[i]),
                                  predicted);
    graph.add(f);
    // Assert that error is null for expected values
    std::vector<Matrix> H(2);
    Vector error = f.unwhitenedError(expected, H);
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(error, gtsam::Vector6::Zero(), 1e-5));
  }

  // Add a prior of 0,0,0 on the first coefficient
  Expression<Coefficient> leaf1(outKeys[0]);
  Expression<ValueType> prior = kindr::minimal::transformationFromHermiteTransformation(leaf1);
  Eigen::Matrix<double, 6, 1> derivatives;
  derivatives << 0,0,0,0,0,0;
  ValueType priorPose(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  graph.add(ExpressionFactor<ValueType>(priorNoiseModel,priorPose,prior));

  // Optimize
  gtsam::LevenbergMarquardtParams params;
  params.setVerbosity("LINEAR");
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials, params).optimize();
  ASSERT_TRUE(expected.equals(result,1e-5));
}

ExpressionFactor<ValueType>
getFactorRelativeMeasurement(const CubicHermiteSE3Curve& curve,
                             Time timeA, Time timeB,
                             ValueType measurement,
                             noiseModel::Diagonal::shared_ptr noiseModel) {
  Expression<ValueType> TA(curve.getValueExpression(timeA));
  Expression<ValueType> TB(curve.getValueExpression(timeB));
  Expression<ValueType> predicted = kindr::minimal::invertAndCompose(TA,TB);
  return ExpressionFactor<ValueType>(noiseModel, measurement, predicted);
}

TEST(CurvesTestSuite, testHermiteSE3RelativeExpression){
  CubicHermiteSE3Curve curve;
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
  curve.extend(times, initials, &outKeys);

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
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, gtsamInitial, params).optimize();

//  gtsamExpected.print("expected:");
//  result.print("result: ");
  ASSERT_TRUE(gtsamExpected.equals(result, 0.5));
}
#endif
