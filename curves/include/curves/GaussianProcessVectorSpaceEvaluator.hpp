#ifndef CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP
#define CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP

#include "Evaluator.hpp"
#include "VectorSpaceConfig.hpp"
#include "GaussianProcessVectorSpaceCurve.hpp"

namespace curves {

class Coefficients;

/// \class GaussianProcessVectorSpaceEvaluator
///
/// Evaluates a GaussianProcessVectorSpaceCurve at a fixed time, for
/// a different set of coefficients.
class GaussianProcessVectorSpaceEvaluator : public Evaluator<VectorSpaceConfig> {

 public:
  /// \brief Parent class
  typedef Evaluator<VectorSpaceConfig> Parent;

  /// \brief The value type of the curve.
  typedef Parent::ValueType ValueType;

  /// \brief The derivative type of the curve.
  typedef Parent::DerivativeType DerivativeType;

  /// \brief Constructor to make a Gaussian process curve evaluator.
  GaussianProcessVectorSpaceEvaluator(const GaussianProcessVectorSpaceCurve& curve, const Time& time);
  virtual ~GaussianProcessVectorSpaceEvaluator();

  /// \brief Get a copy of the vector of keys
  virtual void getKeys(std::vector<Key> *outKeys) const;

  /// \brief Append a copy of the vector of keys to outKeys
  virtual void appendKeys(std::vector<Key> *outKeys) const;

  /// \brief Get an iterator referring to the first key
  virtual std::vector<Key>::const_iterator keyBegin() const;

  /// \brief Get an iterator referring to the past-the-end key in the vector container.
  virtual std::vector<Key>::const_iterator keyEnd() const;

  /// Evaluate the ambient space of the curve (functional form) by specifying new coefficients.
  virtual ValueType evaluate(const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivative(unsigned derivativeOrder,
                                             const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the ambient space of the curve (functional form).
  virtual Eigen::VectorXd evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                              std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the curve derivatives (functional form).
  virtual Eigen::VectorXd evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                        const std::vector<Coefficient>& coefficients,
                                                        std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluate(const Coefficients& coefficients) const;

  /// Get the curve Jacobians.
  /// This is the main interface for GTSAM
  virtual void getJacobians(unsigned derivativeOrder,
                            const Coefficients& coefficients,
                            const Eigen::MatrixXd& chainRule,
                            const std::vector<Eigen::MatrixXd*>& jacobians) const;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluateDerivative(unsigned derivativeOrder,
                                       const Coefficients& coefficients) const;

 protected:
  /// Vector of coefficient keys associated with this evaluator
  std::vector<Key> keys_;

  /// Fixed query time of the evaluation
  Time queryTime_;

  /// Key times associated with the coefficient keys
  std::vector<Time> keyTimes_;

  /// Reference to the curve's prior function
  boost::shared_ptr<const GaussianProcessVectorSpacePrior> prior_;
};

} // namespace curves

#endif /* CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_EVALUATOR_HPP */
