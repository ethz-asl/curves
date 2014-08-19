#include <curves/KeyGenerator.hpp>
#include <boost/thread.hpp>

namespace curves {

size_t KeyGenerator::getNextKey() {
  static size_t key = 0;
  static boost::mutex mutex;
  boost::mutex::scoped_lock lock(mutex);
  return ++key;
}


} // namespace curves
