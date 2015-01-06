/*
 * @file EvaluatorBase.hpp
 * @date Aug 19, 2014
 * @author Paul Furgale
 */

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

  // TODO These imply that derived classes will easily be able to
  // generate vector iterators. Perhaps there is a way to remove this restriction?
  virtual std::vector<Key>::const_iterator keyBegin() const = 0;

  virtual std::vector<Key>::const_iterator keyEnd() const = 0;

  /// \brief Get keys for the coefficients that this evaluator uses.
  ///        This method appends the keys to the vector
  virtual void appendKeys(std::vector<Key> *outKeys) const = 0;

};

} // namespace curves

#endif /* CURVES_EVALUATOR_HPP */
