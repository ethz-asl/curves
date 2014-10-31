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
  typedef typename Parent::EvaluatorType EvaluatorType;
  typedef typename Parent::EvaluatorTypePtr EvaluatorTypePtr;

  VectorSpaceCurve() : dimension_(N) { }
  virtual ~VectorSpaceCurve() {}

  /// \brief Get the dimension of this curve
  virtual size_t dim() const {
    return N;
  }
 
 private:
  /// The dimension of the vector space.
  size_t dimension_;
};

} // namespace curves


#endif /* CT_VECTOR_SPACE_CURVE_HPP */
