/*
 * RosMultiDOFJointTrajectoryInterface.cpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Peter Fankhauser
 */

#include <curves_ros/RosMultiDOFJointTrajectoryInterface.hpp>
#include <ros/ros.h>

// STD
#include <vector>

// Kindr
#include "kindr_ros/kindr_ros.hpp"

namespace curves {

RosMultiDOFJointTrajectoryInterface::RosMultiDOFJointTrajectoryInterface() :
    CubicHermiteSE3Curve()
{

}

RosMultiDOFJointTrajectoryInterface::~RosMultiDOFJointTrajectoryInterface()
{

}

bool RosMultiDOFJointTrajectoryInterface::fromMessage(
    const trajectory_msgs::MultiDOFJointTrajectory& message, const std::string& jointName)
{
  size_t nJoints = message.joint_names.size();
  size_t j = 0;
  for (; j < nJoints; ++j) {
    if (message.joint_names[j] == jointName) break;
    if (j == nJoints - 1) return false; // Joint name not found.
  }

  std::vector<Time> times;
  std::vector<CubicHermiteSE3Curve::ValueType> values;

  for (const auto& point : message.points) {
    times.push_back(ros::Duration(point.time_from_start).toSec());
    CubicHermiteSE3Curve::ValueType pose;
    kindr_ros::convertFromRosGeometryMsg(point.transforms[j], pose);
    values.push_back(pose);
  }

  // TODO Make this work also with velocities and accelerations.

  fitCurve(times, values);

  return true;
}

}
