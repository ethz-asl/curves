/*
 * PolynomialSplineScalarCurve.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Péter Fankhauser, Dario Bellicoso
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once


// stl
#include <string>
#include <vector>

// glog
#include <glog/logging.h>

// curves
#include "curves/Curve.hpp"
#include "curves/ScalarCurveConfig.hpp"
#include "curves/PolynomialSplineContainer.hpp"
#include "curves/polynomial_splines_containers.hpp"

namespace curves {

template<typename SplineContainerType_>
class PolynomialSplineScalarCurve : public Curve<ScalarCurveConfig>
{
 public:
  typedef Curve<ScalarCurveConfig> Parent;
  typedef typename Parent::ValueType ValueType;
  typedef typename Parent::DerivativeType DerivativeType;

  using SplineContainerType = SplineContainerType_;

  PolynomialSplineScalarCurve()
      : Parent(),
        container_(),
        minTime_(0.0)
  {
  }

  virtual ~PolynomialSplineScalarCurve() {

  }

  virtual void print(const std::string& str = "") const
  {
    const double minTime = getMinTime();
    const double maxTime = getMaxTime();
    double timeAtEval = minTime;
    int nPoints = 15;
    double timeDiff = (maxTime-minTime)/(nPoints-1);

    for (int i=0;i<nPoints;i++) {
      double firstDerivative;
      double secondDerivative;
      double value;
      evaluate(value, timeAtEval);
      evaluateDerivative(firstDerivative, timeAtEval, 1);
      evaluateDerivative(secondDerivative, timeAtEval, 2);
      printf("t: %lf, x: %lf dx: %lf dxx: %lf\n",
            timeAtEval,
            value,
            firstDerivative,
            secondDerivative);
      timeAtEval += timeDiff;
    }
  }

  virtual Time getMinTime() const
  {
    return minTime_;
  }

  virtual Time getMaxTime() const
  {
    return container_.getContainerDuration() + minTime_;
  }

  virtual bool evaluate(ValueType& value, Time time) const
  {
    time -= minTime_;
    value = container_.getPositionAtTime(time);
    return true;
  }

  virtual bool evaluateDerivative(DerivativeType& value, Time time, unsigned derivativeOrder) const
  {
    time -= minTime_;
    switch (derivativeOrder) {
      case(1): {
        value = container_.getVelocityAtTime(time);
      } break;

      case(2): {
        value = container_.getAccelerationAtTime(time);
      } break;

      default:
        return false;
    }

    return true;
  }

  virtual void extend(const std::vector<Time>& times, const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys)
  {
    throw std::runtime_error("extend is not yet implemented!");
  }

  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL)
  {
    container_.setData(times, values, 0.0, 0.0, 0.0, 0.0);
    minTime_ = times.front();
  }

  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values,
                        double initialVelocity, double initialAcceleration,
                        double finalVelocity, double finalAcceleration,
                        std::vector<Key>* outKeys = NULL)
  {
    container_.setData(times, values, initialVelocity, initialAcceleration, finalVelocity,
                       finalAcceleration);
    minTime_ = times.front();
  }

  virtual void fitCurve(const std::vector<SplineOptions>& optionList,
                        std::vector<Key>* outKeys = NULL)
  {
    container_.reserveSplines(optionList.size());
    for (const auto& options : optionList) {
      container_.addSpline(PolynomialSplineQuintic(options));
    }
    minTime_ = 0.0;
  }

  virtual void clear()
  {
    container_.reset();
    minTime_ = 0.0;
  }

  virtual void transformCurve(const ValueType T)
  {
    CHECK(false) << "Not implemented";
  }

 protected:
  SplineContainerType container_;
  Time minTime_;
};

using PolynomialSplineQuinticScalarCurve = PolynomialSplineScalarCurve<PolynomialSplineContainerQuintic>;

} /* namespace curves */

