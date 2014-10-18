#ifndef CURVES_GAUSSIAN_PROCESS_PRIOR_FACTOR_EVALUATOR_HPP
#define CURVES_GAUSSIAN_PROCESS_PRIOR_FACTOR_EVALUATOR_HPP

#include "EvaluatorBase.hpp"

namespace curves {

class Coefficients;

class GaussianProcessPriorFactorEvaluator : public EvaluatorBase {

 public:

  GaussianProcessPriorFactorEvaluator();
  virtual ~GaussianProcessPriorFactorEvaluator();

  virtual void getKeys(std::vector<Key>* outKeys) const;

  virtual std::vector<Key>::const_iterator keyBegin() const;

  virtual std::vector<Key>::const_iterator keyEnd() const;

  /// \brief Get keys for the coefficients that this evaluator uses.
  ///        This method appends the keys to the vector
  virtual void appendKeys(std::vector<Key> *outKeys) const;

  /// Evaluate the unwhitened error of prior factor
  virtual Eigen::VectorXd unwhitenedError(const Coefficients& coefficients) const = 0;

  /// Get the prior Jacobians.
  /// This is the main interface for GTSAM
  virtual Eigen::VectorXd unwhitenedErrorAndJacobians(const Coefficients& coefficients,
                                                      const std::vector<Eigen::MatrixXd*>& jacobians) const = 0;

  virtual Eigen::MatrixXd getCovariance() const = 0;
  virtual Eigen::MatrixXd getInformation() const = 0;

 protected:

  std::vector<Key> keys_;

};

} // namespace curves


#endif /* CURVES_GAUSSIAN_PROCESS_PRIOR_FACTOR_EVALUATOR_HPP */
