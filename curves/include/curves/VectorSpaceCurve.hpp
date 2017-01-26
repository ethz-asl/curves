/*
 * ScalarCurveConfig.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Paul Furgale, Renaud Dube, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include "curves/Curve.hpp"
#include "curves/VectorSpaceConfig.hpp"

namespace curves {

template <int N>
class VectorSpaceCurve : public Curve<VectorSpaceConfig<N> >
{
 public:
  typedef Curve<VectorSpaceConfig<N> > Parent;
  typedef typename Parent::ValueType ValueType;
  typedef typename Parent::DerivativeType DerivativeType;

  VectorSpaceCurve() : dimension_(N) { }
  virtual ~VectorSpaceCurve() {}

  /// \brief Get the dimension of this curve
  size_t dim() const {
    return N;
  }

 private:
  /// The dimension of the vector space.
  size_t dimension_;
};

} // namespace
