/*
 * RosJointTrajectoryInterface.cpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Peter Fankhauser
 */

#include <curves_ros/RosJointTrajectoryInterface.hpp>
#include <ros/ros.h>

// STD
#include <vector>

namespace curves {

RosJointTrajectoryInterface::RosJointTrajectoryInterface() :
    PolynomialSplineQuinticScalarCurve()
{

}

RosJointTrajectoryInterface::~RosJointTrajectoryInterface()
{

}

bool RosJointTrajectoryInterface::fromMessage(const trajectory_msgs::JointTrajectory& message,
                                              const std::string& jointName)
{
  size_t nJoints = message.joint_names.size();
  size_t j = 0;
  for (; j < nJoints; ++j) {
    if (message.joint_names[j] == jointName) break;
    if (j == nJoints - 1) return false; // Joint name not found.
  }

  std::vector<Time> times;
  std::vector<PolynomialSplineQuinticScalarCurve::ValueType> values;

  for (const auto& point : message.points) {
    times.push_back(ros::Duration(point.time_from_start).toSec());
    values.push_back(point.positions[j]);
  }

  // TODO Make this work also with velocities and accelerations.

  fitCurve(times, values);

  return true;
}

}
