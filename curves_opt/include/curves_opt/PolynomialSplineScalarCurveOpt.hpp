/*
 * PolynomialSplineScalarCurve.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Péter Fankhauser, Dario Bellicoso
 */

#pragma once

// curves
#include <curves/PolynomialSplineScalarCurve.hpp>
#include <curves/polynomial_splines.hpp>

// curves_opt
#include <curves_opt/polynomial_splines_containers.hpp>

namespace curves {

template<typename SplineType>
class PolynomialSplineScalarCurveOpt : public PolynomialSplineScalarCurve<SplineType> {
 private:
  using Base = PolynomialSplineScalarCurve<SplineType>;

 public:
  using ValueType = typename Base::Parent::ValueType;
  using DerivativeType = typename Base::Parent::DerivativeType;

  PolynomialSplineScalarCurveOpt() : Base() { }
  virtual ~PolynomialSplineScalarCurveOpt() { }

  void fitCurveOptimized(
      const std::vector<Time>& times, const std::vector<ValueType>& values,
      double initialVelocity, double initialAcceleration,
      double finalVelocity, double finalAcceleration,
      double weightMinAccel)
  {
    this->container_.setDataOptimized(
        times, values, initialVelocity, initialAcceleration, finalVelocity,
        finalAcceleration, weightMinAccel);
    this->minTime_ = times.front();
  }

};

using PolynomialSplineQuinticScalarCurveOpt = PolynomialSplineScalarCurveOpt<PolynomialSplineContainerOptQuintic>;

} /* namespace curves */
