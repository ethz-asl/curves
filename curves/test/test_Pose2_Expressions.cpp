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
const double fd_step = 1e-6;

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
    //    std::cout << "pval[2] " << pval[2] << std::endl;
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

  for(int i=0; i < N_TEST_ITERATIONS; ++i) {
    T1val = makeRandomSE2();
    T2val = makeRandomSE2();

    SE2 T2valZeroDiffTheta(T2val.x(),T2val.y(),T1val.theta());

    values.update(1, T1val);
    values.update(2, T2valZeroDiffTheta);
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
    Tval = makeRandomSE2();
    values.update(1, Tval);
    //    std::cout << "Tval.theta() " << Tval.theta() << std::endl;
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(logT, values, fd_step, tolerance);
    //    std::cout << "Completed one testing" << std::endl;

    SE2 Tval_theta_0(Tval.x(), Tval.y(), 2e-16);
    values.update(1, Tval_theta_0);
//    std::cout << "Tval_theta_0 x" << Tval_theta_0.x() << " y " << Tval_theta_0.y() << std::endl;
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

TEST(Pose2ExpressionTestSuite, testSlerpStepByStep) {
  SE2 T0val = makeRandomSE2();
  SE2 T1val = makeRandomSE2();

  // Create some values
  Values values;
  values.insert(1, T0val);
  values.insert(2, T1val);

  ETransformation T0(1);
  ETransformation T1(2);

  std::vector<double> alphas;
  alphas.push_back(1e-5);
  alphas.push_back(1.0 - 1e-5);
  alphas.push_back(0.25);
  alphas.push_back(0.75);

  for(int i=0; i < N_TEST_ITERATIONS; ++i) {
    T0val = makeRandomSE2();
    T1val = makeRandomSE2();

    SE2 T1valZeroDiffTheta(T1val.x(),T1val.y(),T0val.theta());

    values.update(1, T0val);
    values.update(2, T1valZeroDiffTheta);

    for (size_t z = 0; z < alphas.size(); ++z) {
      ETransformation T0_inv = inverse(T0);
      ETransformation T0_inv_T1 = compose(T0_inv, T1);
      Expression<Vector3> log_T0_inv_T1 = transformationLog(T0_inv_T1);
      Expression<Vector3> scaled_log_T0_inv_T1 = vectorScaling(log_T0_inv_T1, alphas[z]);
      ETransformation exp_scaled_log_T0_inv_T1 = transformationExp(scaled_log_T0_inv_T1);
      ETransformation slerp = compose(T0, exp_scaled_log_T0_inv_T1);

//      std::cout << "Testing T0_inv" << std::endl;
      testExpressionJacobians(T0_inv, values, fd_step, tolerance);
//      std::cout << "Testing T0_inv_T1" << std::endl;
      testExpressionJacobians(T0_inv_T1, values, fd_step, tolerance);
//      std::cout << "Testing log_T0_inv_T1" << std::endl;
      testExpressionJacobians(log_T0_inv_T1, values, fd_step, tolerance);
//      std::cout << "Testing scaled_log_T0_inv_T1" << std::endl;
      testExpressionJacobians(scaled_log_T0_inv_T1, values, fd_step, tolerance);
//      std::cout << "Testing exp_scaled_log_T0_inv_T1" << std::endl;
      testExpressionJacobians(exp_scaled_log_T0_inv_T1, values, fd_step, tolerance);
//      std::cout << "Testing slerp" << std::endl;
      testExpressionJacobians(slerp, values, fd_step, tolerance);

    }
  }
}

TEST(Pose2ExpressionTestSuite, testTransformFromCurve) {
  SE2 Tval = makeRandomSE2();

  Vector2 v;
  v.setRandom();
  Point2 Pval(v);

  // Create some values
  Values values;
  values.insert(1, Tval);
  values.insert(2, Pval);

  ETransformation T(1);
  Expression<Point2> P(2);

  Expression<Point2> PinW = transformFromCurve(T, P);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    values.update(1, makeRandomSE2());
    v.setRandom();
    values.update(2, Point2(v));

    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(PinW, values, fd_step, tolerance);
  }
}

TEST(Pose2ExpressionTestSuite, testDistanceBetweenPoints) {
  Vector2 v;
  v.setRandom();
  Point2 Pval1(v);
  v.setRandom();
  Point2 Pval2(v);

  // Create some values
  Values values;
  values.insert(1, Pval1);
  values.insert(2, Pval2);

  Expression<Point2> P1(1);
  Expression<Point2> P2(2);

  Expression<double> dist = distanceBetweenPoints(P1, P2);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    v.setRandom();
    values.update(1, Point2(v));
    v.setRandom();
    values.update(2, Point2(v));

    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(dist, values, fd_step, tolerance);
  }

  // todo the following test fails which motivated the use of pointsSubstraction
  // instead of distance between points
  //  {
  //    values.update(1, Point2(0,0));
  //    values.update(2, Point2(0,0));
  //    testExpressionJacobians(dist, values, fd_step, tolerance);
  //  }
}

TEST(Pose2ExpressionTestSuite, testPointsSubtraction) {
  Vector2 v;
  v.setRandom();
  Point2 Pval1(v);
  v.setRandom();
  Point2 Pval2(v);

  // Create some values
  Values values;
  values.insert(1, Pval1);
  values.insert(2, Pval2);

  Expression<Point2> P1(1);
  Expression<Point2> P2(2);

  Expression<Point2> sub = pointsSubtraction(P1, P2);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    v.setRandom();
    values.update(1, Point2(v));
    v.setRandom();
    values.update(2, Point2(v));

    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(sub, values, fd_step, tolerance);
  }

  {
    values.update(1, Point2(0,0));
    values.update(2, Point2(0,0));
    testExpressionJacobians(sub, values, fd_step, tolerance);
  }
}

//TEST(Pose2ExpressionTestSuite, testPose2AsPoint2) {
//  SE2 Tval = makeRandomSE2();
//
//  // Create some values
//  Values values;
//  values.insert(1, Tval);
//
//  ETransformation T(1);
//
//  Expression<Point2> point = pose2AsPoint2(T);
//
//  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
//    values.update(1, makeRandomSE2());
//
//    SCOPED_TRACE("Testing Expression Jacobians.");
//    testExpressionJacobians(point, values, fd_step, tolerance);
//  }
//}

TEST(Pose2ExpressionTestSuite, testPoseRange) {
  SE2 Tval1 = makeRandomSE2();
  SE2 Tval2 = makeRandomSE2();
  // Create some values
  Values values;
  values.insert(1, Tval1);
  values.insert(2, Tval2);

  ETransformation T1(1);
  ETransformation T2(2);

  Expression<double> range = poseRange(T1, T2);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    values.update(1, makeRandomSE2());
    values.update(2, makeRandomSE2());
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(range, values, fd_step, tolerance);
  }
}

TEST(Pose2ExpressionTestSuite, testPointRange) {
  SE2 Tval = makeRandomSE2();
  Vector2 v;
  v.setRandom();
  Point2 Pval(v);
  // Create some values
  Values values;
  values.insert(1, Tval);
  values.insert(2, Pval);

  ETransformation T(1);
  Expression<Point2> P(2);

  Expression<double> range = pointRange(T, P);

  for (int i = 0; i < N_TEST_ITERATIONS; ++i) {
    values.update(1, makeRandomSE2());

    v.setRandom();
    values.update(2, Point2(v));
    SCOPED_TRACE("Testing Expression Jacobians.");
    testExpressionJacobians(range, values, fd_step, tolerance);
  }
}


