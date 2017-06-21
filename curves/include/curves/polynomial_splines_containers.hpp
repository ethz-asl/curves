/*
 * polynomial_splines_containers.hpp
 *
 *  Created on: Mar 7, 2017
 *      Author: Dario Bellicoso
 */

#pragma once

// curves
#include "curves/PolynomialSplineContainer.hpp"

namespace curves {

using PolynomialSplineContainerLinear    = PolynomialSplineContainer<1>;
using PolynomialSplineContainerQuadratic = PolynomialSplineContainer<2>;
using PolynomialSplineContainerCubic     = PolynomialSplineContainer<3>;
using PolynomialSplineContainerQuartic   = PolynomialSplineContainer<4>;
using PolynomialSplineContainerQuintic   = PolynomialSplineContainer<5>;

}
