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

#include "test_Helpers.hpp"

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
  testExpressionJacobians(eTransform, values, fd_step, tolerance);
  testExpressionJacobians(eAngVel, values, fd_step, tolerance);
  testExpressionJacobians(eVel, values, fd_step, tolerance);

}

// test for correct keys and evaluation function in Slerp SE3 curves
TEST(CurvesTestSuite, testCubicHermiteEvaluation) {

  // setup three SE3 objects
  const double t[] = {0, 100, 200};
  std::vector<Time> times(t,t+3);
  std::vector<ValueType> values;
  for (size_t i = 0; i<3; ++i) {
    SO3 Rval;
    Eigen::Vector3d tval;
    Rval.setRandom();
    tval.setRandom();
    ValueType T(Rval, tval);
    values.push_back(T);
  }

  // interpolate curve
  CubicHermiteSE3Curve curve;
  curve.fitCurve(times, values);

  // fill retrieved coefficients in gtsam values container
  Values gtsamValues;
  curve.initializeGTSAMValues(&gtsamValues);

  // ensure that getValuexpression and evaluate produce same values
  for (size_t i = 0; i < 200; ++i) {
    std::vector<std::string> pose, pose2;
    SE3 eval = curve.evaluate(i);
    Expression<ValueType> expressionC = curve.getValueExpression(i);
    ValueType resultC = expressionC.value(gtsamValues);

    ValueType negation = eval * resultC.inverted();
    gtsam::Vector4 noRot;
    noRot << 1,0,0,0;

    EXPECT_TRUE(EIGEN_MATRIX_NEAR(negation.getPosition(), gtsam::Vector3::Zero(), 1e-5));
    EXPECT_TRUE(EIGEN_MATRIX_NEAR(negation.getRotation().vector(), noRot, 1e-5));
  }
}

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
  curve.fitCurve(times, coefficients, &outKeys);

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
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials).optimize();
  ASSERT_TRUE(expected.equals(result,1e-5));
}
