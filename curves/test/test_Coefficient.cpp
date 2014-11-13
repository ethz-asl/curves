/*
 * @file test_Coefficient.cpp
 * @date Aug 17, 2014
 * @author Paul Furgale
 */

#include <gtest/gtest.h>
#include <curves/Coefficient.hpp>
#include <curves/VectorSpaceCoefficientImplementation.hpp>

TEST(CurvesTestSuite, testCoefficientDefaultConstructor) {
  using namespace curves;
  Coefficient c1;
  ASSERT_EQ(0, c1.dim());  
  ASSERT_EQ(0, c1.ambientDim());

  Coefficient c2;
  c2 = c1;
  ASSERT_EQ(0, c2.dim());  
  ASSERT_EQ(0, c2.ambientDim());

  Coefficient c3(c1);
  ASSERT_EQ(0, c3.dim());  
  ASSERT_EQ(0, c3.ambientDim());

}

TEST(CurvesTestSuite, testCoefficientCopyConstructor) {
  using namespace curves;
  Eigen::VectorXd v = Eigen::VectorXd::Random(3);

  Coefficient c1(v);
  ASSERT_EQ(3, c1.dim());  
  ASSERT_EQ(3, c1.ambientDim());

  Coefficient c2;
  c2 = c1;
  ASSERT_EQ(3, c2.dim());  
  ASSERT_EQ(3, c2.ambientDim());

  Coefficient c3(c1);
  ASSERT_EQ(3, c3.dim());  
  ASSERT_EQ(3, c3.ambientDim());

  // \todo PTF Replace with Eigen assert macros from https://github.com/ethz-asl/eigen_checks
  for(unsigned i = 0; i < 3; ++i) {
    ASSERT_EQ(v[i], c1[i]) << "Difference at index " << i;
    ASSERT_EQ(v[i], c2[i]) << "Difference at index " << i;
    ASSERT_EQ(v[i], c3[i]) << "Difference at index " << i;
  }

}


TEST(CurvesTestSuite, testCoefficientSetVector) {
  using namespace curves;
  Eigen::VectorXd v1 = Eigen::VectorXd::Random(3);
  Coefficient c1(v1);
    
  // \todo PTF Replace with Eigen assert macros from https://github.com/ethz-asl/eigen_checks
  for(unsigned i = 0; i < 3; ++i) {
    ASSERT_EQ(v1[i], c1[i]) << "Difference at index " << i;
  }

  ASSERT_TRUE(c1.equals(Coefficient(v1)));

  Eigen::VectorXd v2 = Eigen::VectorXd::Random(3);
  c1.setValue(v2);

  // \todo PTF Replace with Eigen assert macros from https://github.com/ethz-asl/eigen_checks
  for(unsigned i = 0; i < 3; ++i) {
    ASSERT_NE(v1[i], c1[i]) << "Difference at index " << i;
    ASSERT_EQ(v2[i], c1[i]) << "Difference at index " << i;
  }

  ASSERT_FALSE(c1.equals(Coefficient(v1)));
  ASSERT_TRUE(c1.equals(Coefficient(v2)));

  Eigen::VectorXd v3 = Eigen::VectorXd::Random(10);
  c1.setVector( v3.head<3>() );
  // \todo PTF Replace with Eigen assert macros from https://github.com/ethz-asl/eigen_checks
  for(unsigned i = 0; i < 3; ++i) {
    ASSERT_NE(v1[i], c1[i]) << "Difference at index " << i;
    ASSERT_NE(v2[i], c1[i]) << "Difference at index " << i;
    ASSERT_EQ(v3[i], c1[i]) << "Difference at index " << i;
  }

  ASSERT_FALSE(c1.equals(Coefficient(v1)));
  ASSERT_FALSE(c1.equals(Coefficient(v2)));
  
}

TEST(CurvesTestSuite, testCoefficientBounds) {
  using namespace curves;
  Eigen::VectorXd v2 = Eigen::VectorXd::Random(2);
  Eigen::VectorXd v3 = Eigen::VectorXd::Random(3);
  Eigen::VectorXd v4 = Eigen::VectorXd::Random(4);

  Coefficient c2(v2);
  Coefficient c3(v3);
  Coefficient c4(v4);

  EXPECT_DEATH(c3.equals(c4),"Check failed");
  EXPECT_DEATH(c4.equals(c3),"Check failed");
  EXPECT_DEATH(c2.equals(c3),"Check failed");
  EXPECT_DEATH(c3.equals(c2),"Check failed");
  EXPECT_DEATH(c2.equals(c4),"Check failed");

  EXPECT_DEATH(c2.retract(v3), "Check failed");
  EXPECT_DEATH(c2.retract(v4), "Check failed");
  EXPECT_DEATH(c2.localCoordinates(c3), "Check failed");
  EXPECT_DEATH(c2.localCoordinates(c4), "Check failed");
  EXPECT_DEATH(c2.setVector(v3), "Check failed");
  EXPECT_DEATH(c2.setVector(v4), "Check failed");
  EXPECT_DEATH(c2.setVector(v4.head<3>()), "Check failed");
  EXPECT_DEATH(c2.setValue(v3), "Check failed");
  EXPECT_DEATH(c2.setValue(v4), "Check failed");

  EXPECT_DEATH(c2[2], "Check failed");
  EXPECT_DEATH(c2[3], "Check failed");
}
