#include <curves/VectorSpaceCurve.hpp>

namespace curves {

VectorSpaceCurve::VectorSpaceCurve(size_t dimension) : 
    dimension_(dimension) { }
VectorSpaceCurve::~VectorSpaceCurve() { }

  /// \brief Get the dimension of this curve
size_t VectorSpaceCurve::dim() const {
  return dimension_;
}


} // namespace curves
