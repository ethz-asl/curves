#include <curves/GaussianProcessPriorFactorEvaluator.hpp>
#include <curves/Coefficients.hpp>

using namespace std;

namespace curves {

  GaussianProcessPriorFactorEvaluator::GaussianProcessPriorFactorEvaluator() {}
  GaussianProcessPriorFactorEvaluator::~GaussianProcessPriorFactorEvaluator() {}

  void GaussianProcessPriorFactorEvaluator::getKeys(std::vector<Key> *outKeys) const {
    CHECK_NOTNULL(outKeys);
    *outKeys = keys_;
  }

  void GaussianProcessPriorFactorEvaluator::appendKeys(std::vector<Key> *outKeys) const {
    CHECK_NOTNULL(outKeys);
    outKeys->insert(outKeys->end(), keys_.begin(), keys_.end());
  }

  std::vector<Key>::const_iterator GaussianProcessPriorFactorEvaluator::keyBegin() const {
    return keys_.begin();
  }

  std::vector<Key>::const_iterator GaussianProcessPriorFactorEvaluator::keyEnd() const {
    return keys_.end();
  }

} // namespace curves
