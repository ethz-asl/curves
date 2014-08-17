#ifndef CT_VECTOR_SPACE_CURVE_HPP
#define CT_VECTOR_SPACE_CURVE_HPP

#include "Curve.hpp"
#include "VectorSpaceEvaluator.hpp"

namespace curves {

class VectorSpaceCurve : public Curve
{
 public:
  VectorSpaceCurve();
  virtual ~VectorSpaceCurve();
  /// \name Methods to grow or shrink the curve.
  ///@{
  
  /// Extend the curve into the future so that it can be evaluated
  /// at the following times. The arguments are smoothing points.
  virtual void extendFront(const std::vector<Time>& times,
                           const std::vector<Eigen::VectorXd>& values) = 0;

  /// Extend the curve into the past so that it can be evaluated
  /// at these times. The times should be less than 
  virtual void extendBack(const std::vector<Time>& times,
                          const std::vector<Eigen::VectorXd>& values) = 0;

  ///@}

  /// \brief Fit the curve to these data points
  ///
  /// Underneath the curve should have some default policy for fitting
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<Eigen::VectorXd>& values) = 0;

  /// \brief Get an evaluator at this time
  virtual VectorSpaceEvaluator::Ptr getEvaluator(Time time) = 0;

};

} // namespace curves


#endif /* CT_VECTOR_SPACE_CURVE_HPP */
