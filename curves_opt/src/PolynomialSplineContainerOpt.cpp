/*
 * PolynomialSplineContainerOpt.cpp
 *
 *  Created on: Sep 2, 2017
 *      Author: Dario Bellicoso
 */

// curves_opt
#include <curves_opt/PolynomialSplineContainerOpt.hpp>

namespace curves {

template<>
bool PolynomialSplineContainerOpt<5>::getAccelerationMinimizerBlock(MinAccMat& mat, const double tf) const {
  const double tf2 = tf*tf;
  const double tf3 = tf2*tf;
  const double tf4 = tf3*tf;
  const double tf5 = tf4*tf;
  const double tf6 = tf5*tf;
  const double tf7 = tf6*tf;

  mat << 400.0/7.0*tf7, 40.0*tf6,       24.0*tf5,       10.0*tf4,
         40.0*tf6,      28.8*tf5,       18.0*tf4,       8.0*tf3,
         24.0*tf5,      18.0*tf4,       12.0*tf3,       6.0*tf2,
         10.0*tf4,      8.0*tf3,        6.0*tf2,        4.0*tf;

  return true;
}

template<>
bool PolynomialSplineContainerOpt<4>::getAccelerationMinimizerBlock(MinAccMat& mat, const double tf) const {
  const double tf2 = tf*tf;
  const double tf3 = tf2*tf;
  const double tf4 = tf3*tf;
  const double tf5 = tf4*tf;

  mat << 28.8*tf5,       18.0*tf4,       8.0*tf3,
         18.0*tf4,       12.0*tf3,       6.0*tf2,
         8.0*tf3,        6.0*tf2,        4.0*tf;

  return true;
}

template<>
bool PolynomialSplineContainerOpt<3>::getAccelerationMinimizerBlock(MinAccMat& mat, const double tf) const {
  const double tf2 = tf*tf;
  const double tf3 = tf2*tf;

  mat << 12.0*tf3,       6.0*tf2,
         6.0*tf2,        4.0*tf;

  return true;
}

template<>
bool PolynomialSplineContainerOpt<2>::getAccelerationMinimizerBlock(MinAccMat& mat, const double tf) const {
  mat << 4.0*tf;
  return true;
}

template <int splineOrder_>
bool PolynomialSplineContainerOpt<splineOrder_>::getAccelerationMinimizerBlock(MinAccMat& mat, const double tf) const {
  std::cout << "[PolynomialSplineContainerOpt::getAccelerationMinimizerBlock] Function has not been implemented so far. Spline order was: " << splineOrder_ << std::endl;
  return false;
}

}


