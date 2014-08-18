#ifndef CT_VECTOR_SPACE_CURVE_HPP
#define CT_VECTOR_SPACE_CURVE_HPP

#include "Curve.hpp"
#include "TypedCurve.hpp"
#include "VectorSpaceConfig.hpp"

namespace curves {

class VectorSpaceCurve : public TypedCurve<VectorSpaceConfig>
{
 public:
  typedef TypedCurve<VectorSpaceConfig> Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::EvaluatorType EvaluatorType;
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;

  VectorSpaceCurve(size_t dimension);
  virtual ~VectorSpaceCurve();

  /// \brief Get the dimension of this curve
  virtual size_t dim() const;
 
  /// \brief Get an evaluator at this time
  //virtual VectorSpaceEvaluator::Ptr getEvaluator(Time time) = 0;
 private:
  /// The dimension of the vector space.
  size_t dimension_;
};

} // namespace curves


#endif /* CT_VECTOR_SPACE_CURVE_HPP */
