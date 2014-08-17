#include <curves/KeyGenerator.hpp>
#include <mutex>

namespace curves {

size_t KeyGenerator::getNextKey() {
  static size_t key = 0;
  static std::mutex mutex;
  std::unique_lock<std::mutex> lock(mutex);
  return ++key;
}


} // namespace curves
