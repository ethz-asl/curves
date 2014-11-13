/*
 * @file KeyGenerator.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale
 */

#ifndef CT_KEYGENERATOR_HPP
#define CT_KEYGENERATOR_HPP

#include <cstddef>

namespace curves {

class KeyGenerator
{
 public:
  
  static size_t getNextKey();
  
};

} // namespace curves


#endif /* CT_KEYGENERATOR_HPP */
