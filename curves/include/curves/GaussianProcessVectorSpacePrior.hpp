#ifndef CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "VectorSpaceCurve.hpp"

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


};

} // namespace curves

#endif /* CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
