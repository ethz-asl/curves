/*
 * polynomial_splines_containers.hpp
 *
 *  Created on: Mar 7, 2017
 *      Author: Dario Bellicoso
 */

#pragma once

// curves
#include "curves_opt/PolynomialSplineContainerOpt.hpp"

namespace curves {

using PolynomialSplineContainerOptLinear    = PolynomialSplineContainerOpt<1>;
using PolynomialSplineContainerOptQuadratic = PolynomialSplineContainerOpt<2>;
using PolynomialSplineContainerOptCubic     = PolynomialSplineContainerOpt<3>;
using PolynomialSplineContainerOptQuartic   = PolynomialSplineContainerOpt<4>;
using PolynomialSplineContainerOptQuintic   = PolynomialSplineContainerOpt<5>;

}
