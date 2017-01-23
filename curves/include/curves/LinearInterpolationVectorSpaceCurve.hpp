/*

 * @file LinearInterpolationVectorSpaceCurve-Inl.hpp
 * @date Oct 31, 2014
 * @author Renaud Dube


#ifndef CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_CURVE_HPP
#define CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_CURVE_HPP

#include "VectorSpaceCurve.hpp"
#include "LocalSupport2CoefficientManager.hpp"
#include "gtsam/nonlinear/NonlinearFactorGraph.h"

namespace curves {
template<int N>
class LinearInterpolationVectorSpaceCurve : public VectorSpaceCurve<N> {
 public:
  typedef VectorSpaceCurve<N> Parent;
  typedef typename Parent::ValueType ValueType;
  typedef typename Parent::DerivativeType DerivativeType;

  /// \brief Initialize with the dimension of the vector space
  LinearInterpolationVectorSpaceCurve();
  virtual ~LinearInterpolationVectorSpaceCurve();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// The first valid time for the curve.
  virtual Time getMinTime() const;

  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const;

  /// Extend the curve so that it can be evaluated at these times.
  /// Try to make the curve fit to the values.
  /// Underneath the curve should have some default policy for fitting.
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values,
                      std::vector<Key>* outKeys = NULL);

  /// \brief Fit a new curve to these data points.
  ///
  /// The existing curve will be cleared.
  /// Underneath the curve should have some default policy for fitting.
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL);

  /// Evaluate the ambient space of the curve.
  virtual ValueType evaluate(Time time) const;

  /// Evaluate the curve derivatives.
  /// linear 1st derivative has following behaviour:
  /// - time is out of bound --> error
  /// - time is between 2 coefficients --> take slope between the 2 coefficients
  /// - time is on coefficient (not last coefficient) --> take slope between coefficient and next coefficients
  /// - time is on last coefficient --> take slope between last-1 and last coefficient
  /// derivatives of order >1 equal 0
  virtual DerivativeType evaluateDerivative(Time time, unsigned derivativeOrder) const;

  /// \brief Get an Expression evaluating the curve at this time
  virtual gtsam::Expression<ValueType> getValueExpression(const Time& time) const;

  virtual gtsam::Expression<DerivativeType> getDerivativeExpression(const Time& time, unsigned derivativeOrder) const;

  /// Initialize a GTSAM values structure with the desired keys
  virtual void initializeGTSAMValues(gtsam::KeySet keys, gtsam::Values* values) const;

  /// Initialize a GTSAM values structure for all keys
  virtual void initializeGTSAMValues(gtsam::Values* values) const;

  // updates the relevant curve coefficients from the GTSAM values structure
  virtual void updateFromGTSAMValues(const gtsam::Values& values);

  virtual void clear();

  void addPriorFactors(gtsam::NonlinearFactorGraph* graph, Time priorTime) const;

  /// \brief Perform a rigid transformation on the left side of the curve
  void transformCurve(const ValueType T);

  virtual Time getTimeAtKey(gtsam::Key key) const;

 private:
  typedef Eigen::Matrix<double,N,1> Coefficient;
  typedef typename LocalSupport2CoefficientManager<Coefficient>::TimeToKeyCoefficientMap TimeToKeyCoefficientMap;
  typedef typename LocalSupport2CoefficientManager<Coefficient>::CoefficientIter CoefficientIter;
  LocalSupport2CoefficientManager<Coefficient> manager_;
};

} // namespace curves

#include "LinearInterpolationVectorSpaceCurve-inl.hpp"

#endif  CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_CURVE_HPP
*/
