/*
 * @file test_Pose2_Expressions.cpp
 * @date Nov 25, 2015
 * @author Renaud Dubé, Abel Gawel
 */

#include <gtest/gtest.h>
#include <curves/Pose2_Expressions.hpp>
#include <kindr/minimal/testing-gtsam.h>
#include <boost/assign/list_of.hpp>

#define DIM 3
#define N_TEST_ITERATIONS 10000

using namespace gtsam;
using namespace std;

const double tolerance = 1e-3;
const double fd_step = 1e-9;


SE2 makeRandomSE2() {
  Vector2 pval;
  pval.setRandom();

  return SE2(((double)rand() / RAND_MAX)*6.2832, Point2(pval));
}

TEST(Pose2ExpressionTestSuite, testExp) {
  Vector3 pval;
  pval.setRandom();

  // Create some values
  Values values;
  values.insert(1, pval);

  EVector3 p(1);
  ETransformation T = transformationExp(p);

  for(int i = 0; i < N_TEST_ITERATIONS; ++i) {
    pval.setRandom();
    values.update(1, pval);
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(T, values, fd_step, tolerance);
  }
}

TEST(Pose2ExpressionTestSuite, testInverse) {
  SE2 Tval = makeRandomSE2();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);
  ETransformation invT = inverse(T);
  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    values.update(1, makeRandomSE2());
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(invT, values, fd_step, tolerance);
  }
}

TEST(Pose2ExpressionTestSuite, testCompose) {
  SE2 T1val = makeRandomSE2();
  SE2 T2val = makeRandomSE2();

  // Create some values
  Values values;
  values.insert(1, T1val);
  values.insert(2, T2val);

  ETransformation T1(1);
  ETransformation T2(2);
  ETransformation T1T2 = compose(T1,T2);

  for(int i=0; i < N_TEST_ITERATIONS; ++i) {
    values.update(1, makeRandomSE2());
    values.update(2, makeRandomSE2());
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(T1T2, values, fd_step, tolerance);
  }
}

TEST(Pose2ExpressionTestSuite, testLog) {
  SE2 Tval = makeRandomSE2();

  // Create some values
  Values values;
  values.insert(1, Tval);

  ETransformation T(1);
  Expression<Vector3> logT = transformationLog(T);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    values.update(1, makeRandomSE2());
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(logT, values, fd_step, tolerance);
  }
}

TEST(Pose2ExpressionTestSuite, testSlerp) {
  SE2 T1val = makeRandomSE2();
  SE2 T2val = makeRandomSE2();

  // Create some values
  Values values;
  values.insert(1, T1val);
  values.insert(2, T2val);

  ETransformation T1(1);
  ETransformation T2(2);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    values.update(1, makeRandomSE2());
    values.update(2, makeRandomSE2());
    {
      ETransformation slerpT1 = slerp(T1, T2, 1e-5);
      SCOPED_TRACE("Testing Expression Jacobians.");
      testExpressionJacobians(slerpT1, values, fd_step, tolerance);
    }
    {
      ETransformation slerpT2 = slerp(T1, T2, 1.0 - 1e-5);
      SCOPED_TRACE("Testing Expression Jacobians.");
      testExpressionJacobians(slerpT2, values, fd_step, tolerance);
    }
    {
      ETransformation slerpTa = slerp(T1, T2, 0.25);
      SCOPED_TRACE("Testing Expression Jacobians.");
      testExpressionJacobians(slerpTa, values, fd_step, tolerance);
    }

    {
      ETransformation slerpTb = slerp(T1, T2, 0.75);
      SCOPED_TRACE("Testing Expression Jacobians.");
      testExpressionJacobians(slerpTb, values, fd_step, tolerance);
    }
  }

}

