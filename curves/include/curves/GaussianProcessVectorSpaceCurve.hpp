#ifndef CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP
#define CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP

#include "VectorSpaceCurve.hpp"
#include "HermiteCoefficientManager.hpp"
#include <curves/GaussianProcessVectorSpacePrior.hpp>

class GaussianProcessVectorSpaceEvaluator; // Forward declaration

namespace curves {

class GaussianProcessVectorSpaceCurve : public VectorSpaceCurve {
 public:
  typedef VectorSpaceCurve::ValueType ValueType;
  typedef VectorSpaceCurve::DerivativeType DerivativeType;
  typedef VectorSpaceCurve::EvaluatorType EvaluatorType;
  typedef VectorSpaceCurve::EvaluatorTypePtr EvaluatorTypePtr;

  /// \brief Initialize with the dimension of the vector space
  GaussianProcessVectorSpaceCurve(boost::shared_ptr<GaussianProcessVectorSpacePrior> prior);
  virtual ~GaussianProcessVectorSpaceCurve();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// \brief Get the coefficients that are active at a certain time.
  /// \todo this should be removed in lieu of something that returns a vector
  virtual void appendCoefficientsAt(const Time& time, std::vector<KeyCoefficientTime*>* outCoefficients) const;

  /// \brief Get the coefficients that are active at a certain time.
  virtual void getCoefficientsAt(const Time& time,
                                 Coefficient::Map* outCoefficients) const;

  /// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
  virtual void getCoefficientsInRange(Time startTime, 
                                      Time endTime, 
                                      Coefficient::Map* outCoefficients) const;

  /// \brief Get all of the curve's coefficients.
  virtual void getCoefficients(Coefficient::Map* outCoefficients) const;

  /// \brief Set a coefficient.
  virtual void setCoefficient(Key key, const Coefficient& value);

  /// \brief Set coefficients.
  virtual void setCoefficients(const Coefficient::Map& coefficients);


  /// The first valid time for the curve.
  virtual Time getMinTime() const;

  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const;

  /// Extend the curve so that it can be evaluated at these times.
  /// Try to make the curve fit to the values.
  /// Underneath the curve should have some default policy for fitting.
  virtual void extend(const std::vector<Time>& times,
                      const std::vector<ValueType>& values);

  /// \brief Fit a new curve to these data points.
  ///
  /// The existing curve will be cleared.
  /// Underneath the curve should have some default policy for fitting.
  virtual void fitCurve(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL);

  /// Evaluate the ambient space of the curve.
  virtual Eigen::VectorXd evaluate(Time time) const;

  /// Evaluate the curve derivatives.
  /// linear 1st derivative has following behaviour:
  /// - time is out of bound --> error
  /// - time is between 2 coefficients --> take slope between the 2 coefficients
  /// - time is on coefficient (not last coefficient) --> take slope between coefficient and next coefficients
  /// - time is on last coefficient --> take slope between last-1 and last coefficient
  /// derivatives of order >1 equal 0
  virtual Eigen::VectorXd evaluateDerivative(Time time, unsigned derivativeOrder) const;

  /// \brief Get an evaluator at this time
  EvaluatorTypePtr getEvaluator(const Time& time) const;

  virtual void setTimeRange(Time minTime, Time maxTime);

  boost::shared_ptr<GaussianProcessVectorSpacePrior> getPrior() const { return prior_; }

 private:
  /// Shared pointer to the prior, which determines the type of iterpolation and available sparsity
  boost::shared_ptr<GaussianProcessVectorSpacePrior> prior_;

  /// The manager type is determined by the prior, as the prior determines the local support of the curve
  /// \todo change this to be a shared pointer to a base class of manager... it is created based on the typdef in the prior
  HermiteCoefficientManager manager_;

};

} // namespace curves

#include "GaussianProcessVectorSpaceEvaluator.hpp"

#endif /* CURVES_GAUSSIAN_PROCESS_VECTOR_SPACE_CURVE_HPP */
