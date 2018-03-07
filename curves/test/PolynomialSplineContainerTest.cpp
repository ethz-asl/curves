/*
 * PolynomialSplineContainerTest.cpp
 *
 *  Created on: Feb 3, 2015
 *      Author: Péter Fankhauser, Christian Gehring
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

// gtest
#include <gtest/gtest.h>

// curves
#include "curves/polynomial_splines_containers.hpp"

TEST(PolynomialSplineContainer, getActiveSplineIndexAtTime)
{
  std::vector<double> knotPos;
  std::vector<double> knotVal;

  knotPos.push_back(0.0);
  knotPos.push_back(1.0);
  knotPos.push_back(2.0);

  knotVal.push_back(0.0);
  knotVal.push_back(1.0);
  knotVal.push_back(2.0);

  curves::PolynomialSplineContainerQuintic polyContainer;
  polyContainer.setData(knotPos, knotVal, 0.0, 0.0, 0.0, 0.0);

  double timeOffset = 0.0;

  EXPECT_EQ(0, polyContainer.getActiveSplineIndexAtTime(-1.0, timeOffset));
  EXPECT_EQ(0, timeOffset);
  EXPECT_EQ(0, polyContainer.getActiveSplineIndexAtTime(0.0, timeOffset));
  EXPECT_EQ(0, timeOffset);
  EXPECT_EQ(0, polyContainer.getActiveSplineIndexAtTime(0.5, timeOffset));
  EXPECT_EQ(0, timeOffset);
  EXPECT_EQ(1, polyContainer.getActiveSplineIndexAtTime(1.0, timeOffset));
  EXPECT_EQ(1.0, timeOffset);
  EXPECT_EQ(1, polyContainer.getActiveSplineIndexAtTime(1.5, timeOffset));
  EXPECT_EQ(1.0, timeOffset);
  EXPECT_EQ(1, polyContainer.getActiveSplineIndexAtTime(2.0, timeOffset));
  EXPECT_EQ(1.0, timeOffset);
  EXPECT_EQ(1, polyContainer.getActiveSplineIndexAtTime(2.5, timeOffset));
  EXPECT_EQ(1.0, timeOffset);
}


TEST(PolynomialSplineContainer, eval) {
  std::vector<double> knotPos;
  std::vector<double> knotVal;

  knotPos.push_back(0.0);
  knotPos.push_back(1.0);
  knotPos.push_back(2.0);

  knotVal.push_back(0.0);
  knotVal.push_back(1.0);
  knotVal.push_back(2.0);

  double initialVelocity = 0.1;
  double initialAcceleration = 0.2;

  double finalVelocity = 0.3;
  double finalAcceleration = 0.4;

  curves::PolynomialSplineContainerQuintic polyContainer;
  polyContainer.setData(knotPos, knotVal, initialVelocity, initialAcceleration, finalVelocity, finalAcceleration);

  for (int i=0; i<knotVal.size()-1; i++) {
    EXPECT_NEAR(knotVal[i], polyContainer.getSpline(i)->getPositionAtTime(0.0), 1e-2 ) << " knot:"  << i;
    EXPECT_NEAR(knotVal[i+1], polyContainer.getSpline(i)->getPositionAtTime(knotPos[i+1]-knotPos[i]), 1e-2) << " knot:"  << i;
  }

    EXPECT_NEAR(knotVal[0], polyContainer.getPositionAtTime(knotPos[0]), 1e-2 );
//  for (int i = 0; i < knotVal.size(); ++i) {
//    EXPECT_NEAR(knotVal[i], polyContainer.getPositionAtTime(knotPos[i]), 1e-2 );
//  }

//  EXPECT_NEAR(initialVelocity, polyContainer.getSpline(0)->getVelocityAtTime(0.0), 1e-2 );
//  EXPECT_NEAR(finalVelocity, polyContainer.getSpline(knotVal.size()-2)->getVelocityAtTime(knotPos[knotPos.size()-1]-knotPos[knotVal.size()-2]), 1e-2 );
//
//  EXPECT_NEAR(initialAcceleration, polyContainer.getSpline(0)->getAccelerationAtTime(0.0), 1e-2 );
//  EXPECT_NEAR(finalAcceleration, polyContainer.getSpline(knotVal.size()-2)->getAccelerationAtTime(knotPos[knotPos.size()-1]-knotPos[knotVal.size()-2]), 1e-2 );

}

TEST(PolynomialSplineContainer, getState)
{
  std::vector<double> knotPos;
  std::vector<double> knotVal;

  knotPos.push_back(0.0);
  knotPos.push_back(1.0);
  knotPos.push_back(2.0);

  knotVal.push_back(0.0);
  knotVal.push_back(1.0);
  knotVal.push_back(2.0);

  curves::PolynomialSplineContainerQuintic polyContainer;
  polyContainer.setData(knotPos, knotVal, 0.0, 0.0, 0.0, 0.0);

  constexpr double containerTime = 0.4;
  polyContainer.setContainerTime(containerTime);

  double timeOffset = 0.0;

  EXPECT_EQ(polyContainer.getPosition(), polyContainer.getPositionAtTime(containerTime));
  EXPECT_EQ(polyContainer.getVelocity(), polyContainer.getVelocityAtTime(containerTime));
  EXPECT_EQ(polyContainer.getAcceleration(), polyContainer.getAccelerationAtTime(containerTime));
}
