#ifndef CURVES_EVALUATOR_HPP
#define CURVES_EVALUATOR_HPP

#include "Coefficient.hpp"
#include "types.hpp"

namespace curves {

class EvaluatorBase {
 public:
  EvaluatorBase();
  virtual ~EvaluatorBase();

  /// \brief Get keys for the coefficients that this evaluator uses.
  ///
  /// This function implies an ordering to the keys. For evaluators
  /// defined downstream, they provide their Jacobians in this order.
  /// This may be revisited as we work out the connection to GTSAM
  virtual void getKeys(std::vector<Key>* outKeys) const = 0;
  
  /// \brief Get the coefficients used by this evaluator.
  virtual void getCoefficients(std::vector<Coefficient>* outCoefficients) const = 0;

};

} // namespace curves

#endif /* CURVES_EVALUATOR_HPP */
