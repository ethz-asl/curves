/*
 * polynomial_splines.hpp
 *
 *  Created on: Mar 7, 2017
 *      Author: Dario Bellicoso
 */

#pragma once

// curves
#include "curves/PolynomialSpline.hpp"

namespace curves {

using PolynomialSplineCubic   = PolynomialSpline<3>;
using PolynomialSplineQuintic = PolynomialSpline<5>;

}
