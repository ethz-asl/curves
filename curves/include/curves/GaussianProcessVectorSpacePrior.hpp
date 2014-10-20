#ifndef CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "VectorSpaceCurve.hpp"

namespace curves {

/// \class GaussianProcessVectorSpacePrior
///
/// Extendable interface for a vector-space, Gaussian-process prior.
class GaussianProcessVectorSpacePrior : public VectorSpaceCurve
{
 public:
  /// \brief Parent class
  typedef VectorSpaceCurve Parent;

  /// \brief The value type of the curve.
  typedef Parent::ValueType ValueType;

  /// \brief The derivative type of the curve.
  typedef Parent::DerivativeType DerivativeType;

  /// \brief The evaluator type of the curve.
  typedef Parent::EvaluatorType EvaluatorType;

  /// \brief The evaluator type pointer.
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;

  /// \todo create base N-local-support coefficient manager
  /// \todo add a CoefficientManager typedef to help GP curves initialize properly

  /// \brief General constructor using the curve dimension
  GaussianProcessVectorSpacePrior(size_t dimension) : VectorSpaceCurve(dimension) { }
  virtual ~GaussianProcessVectorSpacePrior() { }

  /// Evaluate the prior at the query time and the key times (associated
  /// with local support), and evaluate the interpolation matrix
  /// associated with the key times (non-zero blocks of K(t)K^{-1})
  virtual Eigen::VectorXd evaluateAndInterpMatrices(Time time, const std::vector<Time>& keyTimes,
                                                    std::vector<Eigen::VectorXd>* outEvalAtKeyTimes,
                                                    std::vector<Eigen::MatrixXd>* outInterpMatrices) const = 0;

  /// Evaluate the derivative of the prior at the query time and the
  /// key times (associated with local support), and evaluate the
  /// interpolation matrix associated with the derivative order and
  /// key times (non-zero blocks of K(t)K^{-1}).
  virtual Eigen::VectorXd evaluateDerivativeAndInterpMatrices(Time time, unsigned derivativeOrder, const std::vector<Time>& keyTimes,
                                                              std::vector<Eigen::VectorXd>* outEvalAtKeyTimes,
                                                              std::vector<Eigen::MatrixXd>* outInterpMatrices) const = 0;

  /// Get the number of key times
  virtual unsigned getNumKeyTimes() const = 0;

 private:

  /// Add a keytime to the prior
  virtual void addKeyTime(const Time& time) = 0;

  /// Add a keytimes to the prior
  virtual void addKeyTimes(const std::vector<Time>& times) = 0;

  /// Clear the keytimes in the prior
  virtual void clearKeyTimes() = 0;

  /// The GP Curve is a friend class so that it can manipulate
  /// the keytimes to mirror it's own internal coefficient times
  friend class GaussianProcessVectorSpaceCurve;
};

} // namespace curves

#endif /* CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
