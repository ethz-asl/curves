#include <curves/SE3Curve.hpp>

namespace curves {

SE3Curve::SE3Curve(){}

SE3Curve::~SE3Curve(){}

size_t SE3Curve::dim() const {
return 6u;
}

}  // namespace curves
