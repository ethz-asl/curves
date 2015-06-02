/*
 * RosMultiDOFJointTrajectoryTranslationInterface.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Peter Fankhauser
 */

#pragma once

// Curves
#include <curves/PolynomialSplineVectorSpaceCurve.hpp>

// ROS
#include <trajectory_msgs/MultiDOFJointTrajectory.h>

// STD
#include <string>

namespace curves {

class RosMultiDOFJointTrajectoryTranslationInterface : public PolynomialSplineQuinticVector3Curve
{
 public:
  RosMultiDOFJointTrajectoryTranslationInterface();
  virtual ~RosMultiDOFJointTrajectoryTranslationInterface();

  /*!
   * Populate spline from a ROS multi DOF joint trajectory message.
   * Considers only the translation.
   * @param message the ROS multi DOF trajectory message.
   * @param jointName the name of the joint to be copied.
   */
  virtual bool fromMessage(const trajectory_msgs::MultiDOFJointTrajectory& message,
                           const std::string& jointName);
};

}  // namespace
