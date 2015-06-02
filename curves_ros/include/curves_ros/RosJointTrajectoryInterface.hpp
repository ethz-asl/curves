/*
 * RosJointTrajectoryInterface.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Peter Fankhauser
 */

#pragma once

// Curves
#include <curves/PolynomialSplineScalarCurve.hpp>

// ROS
#include <trajectory_msgs/JointTrajectory.h>

// STD
#include <string>

namespace curves {

class RosJointTrajectoryInterface : public PolynomialSplineQuinticScalarCurve
{
 public:
  RosJointTrajectoryInterface();
  virtual ~RosJointTrajectoryInterface();

  /*!
   * Populate spline from a ROS joint trajectory message.
   * @param message the ROS trajectory message.
   * @param jointName the name of the joint to be copied.
   */
  bool fromMessage(const trajectory_msgs::JointTrajectory& message, const std::string& jointName);
};

}  // namespace
