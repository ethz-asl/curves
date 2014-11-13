/*
 * @file KeyGenerator.cpp
 * @date Aug 19, 2014
 * @author Paul Furgale
 */

#include <curves/KeyGenerator.hpp>
#include <boost/thread.hpp>

namespace curves {

size_t KeyGenerator::getNextKey() {
  static size_t key = 0;
  static boost::mutex mutex;
  boost::lock_guard<boost::mutex> guard(mutex);
  return ++key;
}


} // namespace curves
