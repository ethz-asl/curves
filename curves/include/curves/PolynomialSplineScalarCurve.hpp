/*
 * PolynomialSplineScalarCurve.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Péter Fankhauser, Dario Bellicoso
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include <string>
#include <vector>

#include "curves/Curve.hpp"
#include "curves/ScalarCurveConfig.hpp"

// Robot Utils
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineContainer.hpp"
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineBase.hpp"
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineQuintic.hpp"
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineCubic.hpp"
#include "robotUtils/timers/ChronoTimer.hpp"

using namespace robotUtils;

namespace curves {

template<typename SplineType>
class PolynomialSplineScalarCurve : public Curve<ScalarCurveConfig>
{
 public:
  typedef Curve<ScalarCurveConfig> Parent;
  typedef typename Parent::ValueType ValueType;
  typedef typename Parent::DerivativeType DerivativeType;

  PolynomialSplineScalarCurve()
      : Curve<ScalarCurveConfig>()
  {

  }

  virtual ~PolynomialSplineScalarCurve()
  {
  }

  virtual void print(const std::string& str = "") const
  {
  }

  virtual Time getMinTime() const
  {
  }

  virtual Time getMaxTime() const
  {
    return container_.getContainerDuration();
  }

  virtual ValueType evaluate(Time time) const
  {
    return container_.getPositionAtTime(time);
  }

  virtual ValueType evaluate(Time time, robotUtils::HighResolutionClockTimer& timer) const
  {
    return container_.getPositionAtTime(time, timer);
  }

  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) const
  {
    switch (derivativeOrder) {
      case(1): {
        return container_.getVelocityAtTime(time);
      } break;

      case(2): {
        return container_.getAccelerationAtTime(time);
      } break;
    }
  }

  virtual void extend(const std::vector<Time>& times, const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys)
  {

  }

  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL)
  {
    container_.setData(times, values, 0.0, 0.0, 0.0, 0.0);
  }

  virtual void fitCurve(const std::vector<PolynomialSplineBase::SplineOpts>& values,
                        std::vector<Key>* outKeys = NULL)
  {
    for (auto& value : values) {
      PolynomialSplineQuintic spline;
      spline.evalCoeffs(value);
      container_.addSpline(spline);
    }
  }

 private:
  PolynomialSplineContainer container_;
};

typedef PolynomialSplineScalarCurve<PolynomialSplineQuintic> PolynomialSplineQuinticScalarCurve;
//typedef PolynomialSplineScalarCurve<PolynomialSplineCubic> PolynomialSplineCubicScalarCurve;
//typedef PolynomialSplineScalarCurve<PolynomialSplineLinear> PolynomialSplineLinearScalarCurve;

} // namespace

