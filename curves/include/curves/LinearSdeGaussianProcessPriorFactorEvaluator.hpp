#ifndef CURVES_LINEAR_SDE_GAUSSIAN_PROCESS_PRIOR_FACTOR_EVALUATOR_HPP
#define CURVES_LINEAR_SDE_GAUSSIAN_PROCESS_PRIOR_FACTOR_EVALUATOR_HPP

#include "GaussianProcessPriorFactorEvaluator.hpp"
#include "LinearSdeGaussianProcessVectorSpacePrior.hpp"

namespace curves {

class Coefficients;

class LinearSdeGaussianProcessPriorFactorEvaluator : public GaussianProcessPriorFactorEvaluator {

 public:

  LinearSdeGaussianProcessPriorFactorEvaluator(const Key& key0,
                                               boost::shared_ptr<LinearSdeGaussianProcessVectorSpacePrior::LinearSdeCoefficient> priorCoeff0);
  LinearSdeGaussianProcessPriorFactorEvaluator(const Key& key0, const Key& key1,
                                               boost::shared_ptr<LinearSdeGaussianProcessVectorSpacePrior::LinearSdeCoefficient> priorCoeff1);
  virtual ~LinearSdeGaussianProcessPriorFactorEvaluator();

  /// Evaluate the unwhitened error of prior factor
  virtual Eigen::VectorXd unwhitenedError(const Coefficients& coefficients) const;

  /// Get the prior Jacobians.
  /// This is the main interface for GTSAM
  virtual Eigen::VectorXd unwhitenedErrorAndJacobians(const Coefficients& coefficients,
                                                      const std::vector<Eigen::MatrixXd*>& jacobians) const;

  virtual Eigen::MatrixXd getCovariance() const;
  virtual Eigen::MatrixXd getInformation() const;

 private:

  boost::shared_ptr<LinearSdeGaussianProcessVectorSpacePrior::LinearSdeCoefficient> priorCoeff_;

};

} // namespace curves


#endif /* CURVES_LINEAR_SDE_GAUSSIAN_PROCESS_PRIOR_FACTOR_EVALUATOR_HPP */
