/*
 * @file test_SE3Coefficient.cpp
 * @date Oct 10, 2014
 * @author Mike Bosse
 */

#include <gtest/gtest.h>
#include <curves/Coefficient.hpp>
#include <curves/SE3CoefficientImplementation.hpp>
#include "kindr/minimal/quat-transformation.h"

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef Eigen::Matrix4d Matrix4d;


TEST(CurvesTestSuite, testSE3CoefficientDefaultConstructor) {
  using namespace curves;
  CoefficientImplementation::Ptr impl(new SE3CoefficientImplementation);
  Coefficient c1(impl);
  ASSERT_EQ(6u, c1.dim());  
  ASSERT_EQ(7u, c1.ambientDim());

  Coefficient c2;
  c2 = c1;
  ASSERT_EQ(6u, c2.dim());  
  ASSERT_EQ(7u, c2.ambientDim());

  Coefficient c3(c1);
  ASSERT_EQ(6u, c3.dim());  
  ASSERT_EQ(7u, c3.ambientDim());

}

TEST(CurvesTestSuite, testSE3CoefficientCopyConstructor) {
  using namespace curves;
  SE3 pose(SO3(SO3::Vector3::Random().eval()),SE3::Position::Random().eval());
  
  Matrix4d mat = pose.getTransformationMatrix();

  CoefficientImplementation::Ptr impl(new SE3CoefficientImplementation);

  Eigen::VectorXd val(7);
  boost::dynamic_pointer_cast<SE3CoefficientImplementation>(impl)->makeValue(mat,&val);
  Coefficient c1(impl,val);

  ASSERT_EQ(6u, c1.dim());  
  ASSERT_EQ(7u, c1.ambientDim());

  Coefficient c2;
  c2 = c1;
  ASSERT_EQ(6u, c2.dim());  
  ASSERT_EQ(7u, c2.ambientDim());

  Coefficient c3(c1);
  ASSERT_EQ(6u, c3.dim());  
  ASSERT_EQ(7u, c3.ambientDim());

  // \todo PTF Replace with Eigen assert macros from https://github.com/ethz-asl/eigen_checks
  for(unsigned i = 0; i < 7u; ++i) {
    ASSERT_EQ(val[i], c1[i]) << "Difference at index " << i;
    ASSERT_EQ(val[i], c2[i]) << "Difference at index " << i;
    ASSERT_EQ(val[i], c3[i]) << "Difference at index " << i;
  }

}


TEST(CurvesTestSuite, testSE3CoefficientRetractAndLocalCoord) {
  using namespace curves;
  for (int iter = 0; iter < 10000; ++iter) {
  SE3 poseA(SO3(SO3::Vector3::Random().eval()),SE3::Position::Random().eval());
  SE3 poseB(SO3(SO3::Vector3::Random().eval()),SE3::Position::Random().eval());

  SE3CoefficientImplementation::Ptr impl(new SE3CoefficientImplementation);
  Eigen::VectorXd val(7);

  boost::dynamic_pointer_cast<SE3CoefficientImplementation>(impl)->makeValue(poseA.getTransformationMatrix(),&val);
  Coefficient Ca(impl,val);

  boost::dynamic_pointer_cast<SE3CoefficientImplementation>(impl)->makeValue(poseB.getTransformationMatrix(),&val);
  Coefficient Cb(impl,val);

  Eigen::VectorXd  delta;
  delta = Ca.localCoordinates(Cb);

  EXPECT_EQ(6u,delta.size()) << "local coordinates should be 6 dof.";

  Coefficient Cc = Ca.retract(delta);
  ASSERT_EQ(6u, Cc.dim()) << "dim should be 6";  
  ASSERT_EQ(7u, Cc.ambientDim()) << "ambient dim should be 7";

  EXPECT_TRUE(Cc.equals(Cb)) << "local coordinates and retract are not compatible at iter " << iter;
  }
}
