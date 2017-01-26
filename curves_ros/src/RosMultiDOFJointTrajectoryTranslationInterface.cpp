/*
 * RosMultiDOFJointTrajectoryTranslationInterface.cpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Peter Fankhauser
 */

#include <curves_ros/RosMultiDOFJointTrajectoryTranslationInterface.hpp>
#include <ros/ros.h>

// STD
#include <vector>

namespace curves {

RosMultiDOFJointTrajectoryTranslationInterface::RosMultiDOFJointTrajectoryTranslationInterface() :
    PolynomialSplineQuinticVector3Curve()
{

}

RosMultiDOFJointTrajectoryTranslationInterface::~RosMultiDOFJointTrajectoryTranslationInterface()
{

}

bool RosMultiDOFJointTrajectoryTranslationInterface::fromMessage(
    const trajectory_msgs::MultiDOFJointTrajectory& message, const std::string& jointName)
{
  size_t nJoints = message.joint_names.size();
  size_t j = 0;
  for (; j < nJoints; ++j) {
    if (message.joint_names[j] == jointName) break;
    if (j == nJoints - 1) return false; // Joint name not found.
  }

  std::vector<Time> times;
  std::vector<PolynomialSplineQuinticVector3Curve::ValueType> values;

  for (const auto& point : message.points) { // TODO Make nicer with ROS to kindr transforms.
    times.push_back(ros::Duration(point.time_from_start).toSec());
    PolynomialSplineQuinticVector3Curve::ValueType position(point.transforms[j].translation.x,
                                                            point.transforms[j].translation.y,
                                                            point.transforms[j].translation.z);
    values.push_back(position);
  }

  // TODO Make this work also with velocities and accelerations.

  fitCurve(times, values);

  return true;
}

}
