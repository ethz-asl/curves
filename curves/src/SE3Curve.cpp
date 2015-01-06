/*
 * @file SE3Curve.cpp
 * @date Oct 03, 2014
 * @author Paul Furgale
 */

#include <curves/SE3Curve.hpp>

namespace curves {

SE3Curve::SE3Curve(){}

SE3Curve::~SE3Curve(){}

size_t SE3Curve::dim() const {
return 6u;
}

}  // namespace curves
