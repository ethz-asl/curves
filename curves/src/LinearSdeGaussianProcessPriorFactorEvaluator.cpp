#include <curves/LinearSdeGaussianProcessPriorFactorEvaluator.hpp>
#include <curves/Coefficients.hpp>
#include <Eigen/LU>

using namespace std;

namespace curves {

  LinearSdeGaussianProcessPriorFactorEvaluator::LinearSdeGaussianProcessPriorFactorEvaluator(const Key& key0,
      boost::shared_ptr<LinearSdeGaussianProcessVectorSpacePrior::LinearSdeCoefficient> priorCoeff0) : priorCoeff_(priorCoeff0) {
    keys_.push_back(key0);
  }
  LinearSdeGaussianProcessPriorFactorEvaluator::LinearSdeGaussianProcessPriorFactorEvaluator(const Key& key0, const Key& key1,
      boost::shared_ptr<LinearSdeGaussianProcessVectorSpacePrior::LinearSdeCoefficient> priorCoeff1) : priorCoeff_(priorCoeff1) {
    keys_.push_back(key0);
    keys_.push_back(key1);
  }

  LinearSdeGaussianProcessPriorFactorEvaluator::~LinearSdeGaussianProcessPriorFactorEvaluator() {}

  Eigen::VectorXd LinearSdeGaussianProcessPriorFactorEvaluator::unwhitenedError(const Coefficients& coefficients) const {
    if (keys_.size() == 1) {
      return priorCoeff_->mean - coefficients.get(keys_[0]).getValue();
    } else { /* keys_.size() == 2 */
      return priorCoeff_->liftedExogenousInput - coefficients.get(keys_[1]).getValue() + priorCoeff_->stateTransitionMatrix * coefficients.get(keys_[0]).getValue();
    }
  }

  Eigen::VectorXd LinearSdeGaussianProcessPriorFactorEvaluator::unwhitenedErrorAndJacobians(const Coefficients& coefficients,
                                                                                            const std::vector<Eigen::MatrixXd*>& jacobians) const {
    CHECK_EQ(jacobians.size(), keys_.size()) << "number of jacobians and keys do not match.";
    CHECK_NOTNULL(jacobians[0]);
    size_t dim_coefficient = coefficients.get(keys_[0]).dim();
    if (keys_.size() == 1) {
      *(jacobians[0]) += Eigen::MatrixXd::Identity(dim_coefficient,dim_coefficient);
      return priorCoeff_->mean - coefficients.get(keys_[0]).getValue();
    } else { /* keys_.size() == 2 */
      CHECK_NOTNULL(jacobians[1]);
      *(jacobians[0]) += priorCoeff_->stateTransitionMatrix;
      *(jacobians[1]) += Eigen::MatrixXd::Identity(dim_coefficient,dim_coefficient);
      return priorCoeff_->liftedExogenousInput - coefficients.get(keys_[1]).getValue() + priorCoeff_->stateTransitionMatrix * coefficients.get(keys_[0]).getValue();
    }
  }

  Eigen::MatrixXd LinearSdeGaussianProcessPriorFactorEvaluator::getCovariance() const {
    return priorCoeff_->inverseLiftedCovarianceMatrix.inverse();
  }

  Eigen::MatrixXd LinearSdeGaussianProcessPriorFactorEvaluator::getInformation() const {
    return priorCoeff_->inverseLiftedCovarianceMatrix;
  }


} // namespace curves
