/*
 * polynomial_splines.hpp
 *
 *  Created on: Mar 7, 2017
 *      Author: dbellicoso
 */

#pragma once

namespace curves {

struct SplineOptions {
  double tf = 0.0;
  double pos0 = 0.0;
  double posT = 0.0;
  double vel0 = 0.0;
  double velT = 0.0;
  double acc0 = 0.0;
  double accT = 0.0;
};

}
