/*
 * RosMultiDOFJointTrajectoryInterface.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Peter Fankhauser
 */

#pragma once

// Curves
#include <curves/CubicHermiteSE3Curve.hpp>

// ROS
#include <trajectory_msgs/MultiDOFJointTrajectory.h>

// STD
#include <string>

namespace curves {

class RosMultiDOFJointTrajectoryInterface : public CubicHermiteSE3Curve
{
 public:
  RosMultiDOFJointTrajectoryInterface();
  virtual ~RosMultiDOFJointTrajectoryInterface();

  /*!
   * Populate spline from a ROS multi DOF joint trajectory message.
   * @param message the ROS multi DOF trajectory message.
   * @param jointName the name of the joint to be copied.
   */
  virtual bool fromMessage(const trajectory_msgs::MultiDOFJointTrajectory& message,
                           const std::string& jointName);
};

}  // namespace
