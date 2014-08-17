#ifndef CURVES_EVALUATOR_HPP
#define CURVES_EVALUATOR_HPP

#include "Coefficient.hpp"
#include "types.hpp"

namespace curves {

class Evaluator {
 public:
  Evaluator();
  virtual ~Evaluator();

  virtual void getKeys(std::vector<Key>& outKeys) = 0;
  
  virtual void getCoefficients(std::vector<Key>& outCoefficients) = 0;

};

} // namespace curves

#endif /* CURVES_EVALUATOR_HPP */
