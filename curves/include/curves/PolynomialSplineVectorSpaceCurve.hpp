/*
 * PolynomialSplineVectorSpaceCurve.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include <string>
#include <vector>
#include <Eigen/Core>
#include <glog/logging.h>

#include "curves/Curve.hpp"
#include "curves/VectorSpaceCurve.hpp"
#include "curves/PolynomialSplineContainer.hpp"
#include "curves/polynomial_splines_containers.hpp"

namespace curves {

template<typename SplineType, int N>
class PolynomialSplineVectorSpaceCurve : public VectorSpaceCurve<N>
{
 public:
  typedef VectorSpaceCurve<N> Parent;
  typedef typename Parent::ValueType ValueType;
  typedef typename Parent::DerivativeType DerivativeType;

  PolynomialSplineVectorSpaceCurve()
      : VectorSpaceCurve<N>(),
        minTime_(0)
  {
    containers_.resize(N);
  }

  virtual ~PolynomialSplineVectorSpaceCurve()
  {
  }

  virtual void print(const std::string& str = "") const
  {
  }

  virtual Time getMinTime() const
  {
    return minTime_;
  }

  virtual Time getMaxTime() const
  {
    return containers_.at(0).getContainerDuration();
  }

  virtual bool evaluate(ValueType& value, Time time) const
  {
    for (size_t i = 0; i < N; ++i) {
      value(i) = containers_.at(i).getPositionAtTime(time);
    }
    return true;
  }

  virtual bool evaluateDerivative(DerivativeType& value, Time time, unsigned derivativeOrder) const
  {

    for (size_t i = 0; i < N; ++i) {
      if (derivativeOrder == 1) {
        value(i) = containers_.at(i).getVelocityAtTime(time);
      }
      else if (derivativeOrder == 2) {
        value(i) = containers_.at(i).getAccelerationAtTime(time);
      }
      else {
        return false;
      }
    }
    return true;
  }

  virtual void extend(const std::vector<Time>& times, const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys)
  {
    throw std::runtime_error("PolynomialSplineVectorSpaceCurve::extend is not yet implemented!");
  }

  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL)
  {
    minTime_ = times.front();
    for (size_t i = 0; i < N; ++i) {
      std::vector<double> scalarValues;
      scalarValues.reserve(times.size());
      for (size_t t = 0; t < times.size(); ++t) scalarValues.push_back(values.at(t)(i));
      containers_.at(i).setData(times, scalarValues, 0.0, 0.0, 0.0, 0.0);
    }
  }

  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        const DerivativeType& initialVelocity,
                        const DerivativeType& initialAcceleration,
                        const DerivativeType& finalVelocity,
                        const DerivativeType& finalAcceleration)
  {
    minTime_ = times.front();
    for (size_t i = 0; i < N; ++i) {
      std::vector<double> scalarValues;
      scalarValues.reserve(times.size());
      for (size_t t = 0; t < times.size(); ++t) scalarValues.push_back(values.at(t)(i));
      containers_.at(i).setData(times, scalarValues, initialVelocity(i), initialAcceleration(i),
                                finalVelocity(i), finalAcceleration(i));
    }
  }

  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values,
                        const std::vector<DerivativeType>& firstDerivatives,
                        const std::vector<DerivativeType>& secondDerivatives,
                        std::vector<Key>* outKeys = NULL)
  {
    minTime_ = times.front();
    for (size_t i = 0; i < N; ++i) {
      std::vector<double> scalarValues, scalarFirstDerivates, scalarSecondDerivates;
      scalarValues.reserve(times.size());
      scalarFirstDerivates.reserve(times.size());
      scalarSecondDerivates.reserve(times.size());
      for (size_t t = 0; t < times.size(); ++t) {
        scalarValues.push_back(values.at(t)(i));
        scalarFirstDerivates.push_back(firstDerivatives.at(t)(i));
        scalarSecondDerivates.push_back(secondDerivatives.at(t)(i));
      }
      // TODO Copy all derivates, right now only first and last are supported.
      containers_.at(i).setData(times, scalarValues,
                                *(scalarFirstDerivates.begin()), *(scalarSecondDerivates.begin()),
                                *(scalarFirstDerivates.end() - 1), *(scalarSecondDerivates.end() - 1));
    }
  }


  virtual void fitCurve(const std::vector<SplineOptions>& values,
                        std::vector<Key>* outKeys = NULL)
  {
    // TODO
    throw std::runtime_error("PolynomialSplineVectorSpaceCurve::fitCurve is not yet implemented!");
  }

  virtual void clear()
  {
    for (size_t i = 0; i < N; ++i) {
      containers_.at(i).reset();
    }
  }

  virtual void transformCurve(const ValueType T)
  {
    CHECK(false) << "Not implemented";
  }

 private:
  std::vector<PolynomialSplineContainerQuintic> containers_;
  Time minTime_;
};

typedef PolynomialSplineVectorSpaceCurve<PolynomialSplineQuintic, 3> PolynomialSplineQuinticVector3Curve;

} // namespace

