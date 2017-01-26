/*
 * KeyGenerator.hpp
 *
 *  Created on: Aug 17, 2014
 *      Author: Paul Furgale, Péter Fankhauser
 *   Institute: ETH Zurich, Autonomous Systems Lab
 */

#pragma once

#include <cstddef>

namespace curves {

class KeyGenerator
{
 public:

  static size_t getNextKey();

};

} // namespace curves
