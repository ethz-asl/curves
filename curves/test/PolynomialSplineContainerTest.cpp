/*****************************************************************************************
* Software License Agreement (BSD License)
*
* Copyright (c) 2014, Christian Gehring, Péter Fankhauser, C. Dario Bellicoso, Stelian Coros
* All rights reserved.
*
*  Redistribution and use in source and binary forms, with or without
*  modification, are permitted provided that the following conditions
*  are met:
*
*   * Redistributions of source code must retain the above copyright
*     notice, this list of conditions and the following disclaimer.
*   * Redistributions in binary form must reproduce the above
*     copyright notice, this list of conditions and the following
*     disclaimer in the documentation and/or other materials provided
*     with the distribution.
*   * Neither the name of Autonomous Systems Lab nor ETH Zurich
*     nor the names of its contributors may be used to endorse or
*     promote products derived from this software without specific
*     prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
*  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
*  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
*  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
*  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
*  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
*  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
*  POSSIBILITY OF SUCH DAMAGE.
*/
/*!
* @file     ContactForceDistributionTest.cpp
* @author   Christian Gehring, Peter Fankhauser
* @date     Feb 3, 2015
* @brief
*/


#include <gtest/gtest.h>

#include "curves/PolynomialSplineContainer.hpp"

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

  curves::PolynomialSplineContainer polyContainer;
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

  curves::PolynomialSplineContainer polyContainer;
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
