/*
 * @file VectorSpaceCurve.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale, Renaud Dube
 */

#ifndef CT_VECTOR_SPACE_CURVE_HPP
#define CT_VECTOR_SPACE_CURVE_HPP

#include "Curve.hpp"
#include "VectorSpaceConfig.hpp"

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

} // namespace curves


#endif /* CT_VECTOR_SPACE_CURVE_HPP */
