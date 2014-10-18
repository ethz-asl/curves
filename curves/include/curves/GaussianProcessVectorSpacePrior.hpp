#ifndef CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "VectorSpaceCurve.hpp"
#include "GaussianProcessPriorFactorEvaluator.hpp"

namespace curves {

class GaussianProcessVectorSpacePrior : public VectorSpaceCurve
{
 public:
  typedef VectorSpaceCurve Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::EvaluatorType EvaluatorType;
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;
  /// \todo create base N-local-support coefficient manager
  // typedef CoefficientManager CurveCoefficientManagerType;

  GaussianProcessVectorSpacePrior(size_t dimension) : VectorSpaceCurve(dimension) { }
  virtual ~GaussianProcessVectorSpacePrior() { }

  /// \name Methods to evaluate the curve
  ///@{

  /// Evaluate the prior at the query time, the key times (associated with local support), and evaluate the interpolation matrix associated with the key times (and belonging to K(t)K^{-1}).
  virtual Eigen::VectorXd evaluateAndInterpMatrices(Time time, const std::vector<Time>& keyTimes,
                                                    const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                    const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const = 0;

  virtual unsigned getNumKeyTimes() const = 0;

  virtual std::vector<boost::shared_ptr<GaussianProcessPriorFactorEvaluator> > getPriorFactors() const = 0;

 private:

  /// Add a keytime to the prior
  virtual void addKeyTime(const Time& time, const Key& assocKey) = 0;

  /// Add a keytimes to the prior
  virtual void addKeyTimes(const std::vector<Time>& times, const std::vector<Key>& assocKeys) = 0;

  /// Clear the keytimes in the prior
  virtual void clearKeyTimes() = 0;

  /// The GP Curve is a friend class so that it can manipulate the keytimes to mirror it's own internal coefficient times
  friend class GaussianProcessVectorSpaceCurve;

};

} // namespace curves

#endif /* CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
