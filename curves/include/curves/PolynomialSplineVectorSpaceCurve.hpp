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

#include "curves/Curve.hpp"
#include "curves/VectorSpaceCurve.hpp"

// Robot Utils
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineContainer.hpp"
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineBase.hpp"
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineQuintic.hpp"
#include "robotUtils/function_approximators/polynomialSplines/PolynomialSplineCubic.hpp"

namespace curves {

template<typename SplineType, int N>
class PolynomialSplineVectorSpaceCurve : public VectorSpaceCurve<N>
{
 public:
  typedef VectorSpaceCurve<N> Parent;
  typedef typename Parent::ValueType ValueType;
  typedef typename Parent::DerivativeType DerivativeType;

  PolynomialSplineVectorSpaceCurve()
      : VectorSpaceCurve<N>()
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
  }

  virtual Time getMaxTime() const
  {
    return containers_.at(0).getContainerDuration();
  }

  virtual ValueType evaluate(Time time) const
  {
    ValueType value;
    for (size_t i = 0; i < N; ++i) {
      value(i) = containers_.at(i).getPositionAtTime(time);
    }
    return value;
  }

  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) const
  {
    ValueType value;
    for (size_t i = 0; i < N; ++i) {
      if (derivativeOrder == 1) {
        value(i) = containers_.at(i).getVelocityAtTime(time);
      } else if (derivativeOrder == 2) {
        value(i) = containers_.at(i).getAccelerationAtTime(time);
      }

    }
    return value;
  }

  virtual void extend(const std::vector<Time>& times, const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys)
  {

  }

  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL)
  {
    for (size_t i = 0; i < N; ++i) {
      std::vector<double> scalarValues;
      scalarValues.reserve(times.size());
      for (size_t t = 0; t < times.size(); ++t) scalarValues.push_back(values.at(t)(i));
      containers_.at(i).setData(times, scalarValues, 0.0, 0.0, 0.0, 0.0);
    }
  }

  virtual void fitCurve(const std::vector<PolynomialSplineBase::SplineOpts>& values,
                        std::vector<Key>* outKeys = NULL)
  {
    // TODO
  }

  virtual void clearCurve() {
    for (size_t i = 0; i < N; ++i) {
      containers_.at(i).reset();
    }
  }

 private:
  std::vector<PolynomialSplineContainer> containers_;
};

typedef PolynomialSplineVectorSpaceCurve<PolynomialSplineQuintic, 3> PolynomialSplineQuinticVector3Curve;

} // namespace

