/*
 * polynomial_splines.cpp
 *
 *  Created on: Mar 7, 2017
 *      Author: Dario Bellicoso
 */

// curves
#include "curves/polynomial_splines_traits.hpp"

namespace curves {
namespace spline_traits {

template<typename Core_, int SplineOrder_>
constexpr unsigned int spline_rep<Core_, SplineOrder_>::numCoefficients;

using spline3 = spline_rep<double, 3>;
constexpr std::array<double, spline3::numCoefficients> spline3::tauZero;
constexpr std::array<double, spline3::numCoefficients> spline3::dtauZero;
constexpr std::array<double, spline3::numCoefficients> spline3::ddtauZero;

using spline5 = spline_rep<double, 5>;
constexpr std::array<double, spline5::numCoefficients> spline5::tauZero;
constexpr std::array<double, spline5::numCoefficients> spline5::dtauZero;
constexpr std::array<double, spline5::numCoefficients> spline5::ddtauZero;

}
}
