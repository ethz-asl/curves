#ifndef CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_CURVE_HPP
#define CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_CURVE_HPP

#include "VectorSpaceCurve.hpp"
#include "HermiteCoefficientManager.hpp"

namespace curves {

class LinearInterpolationVectorSpaceCurve : public VectorSpaceCurve {
 public:
  typedef VectorSpaceCurve::ValueType ValueType;
  typedef VectorSpaceCurve::DerivativeType DerivativeType;
  typedef VectorSpaceCurve::EvaluatorType EvaluatorType;
  typedef VectorSpaceCurve::EvaluatorTypePtr EvaluatorTypePtr;

  /// \brief Initialize with the dimension of the vector space
  LinearInterpolationVectorSpaceCurve(size_t dimension);
  virtual ~LinearInterpolationVectorSpaceCurve();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;
  
  /// \name Methods to get and set coefficients.
  ///@{

  /// \brief Get the coefficients that are active at a certain time.
  virtual void getCoefficientsAt(Time time, 
                                 Coefficient::Map& outCoefficients) const;

  /// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
  virtual void getCoefficientsInRange(Time startTime, 
                                      Time endTime, 
                                      Coefficient::Map& outCoefficients) const;

  /// \brief Get all of the curve's coefficients.
  virtual void getCoefficients(Coefficient::Map& outCoefficients) const;
  
  /// \brief Set a coefficient.
  virtual void setCoefficient(Key key, const Coefficient& value);

  /// \brief Set coefficients.
  virtual void setCoefficients(Coefficient::Map& coefficients);
  ///@}


  /// \name Methods to deal with time.
  ///@{

  /// The first valid time for the curve.
  virtual Time getBackTime() const;
  
  /// The one past the last valid time for the curve.
  virtual Time getFrontTime() const;

  ///@}

  /// \name Methods to grow or shrink the curve.

  /// Extend the curve into the future so that it can be evaluated
  /// at this time.
  virtual void extendFront(Time time);

  /// Extend the curve into the past so that it can be evaluated
  /// at this time.
  virtual void extendBack(Time time);

  /// Retract the curve in the front. The caller guarantees that
  /// they will never try to evaluate this time or higher again.
  virtual void retractFront(Time time);

  /// Retract the curve in the back. The caller guarantees that
  /// they will never try to evaluate this time or lower again.
  virtual void retractBack(Time time);

  /// Extend the curve into the future so that it can be evaluated
  /// at the following times. The arguments are smoothing points.
  virtual void extendFront(const std::vector<Time>& times,
                           const std::vector<Eigen::VectorXd>& values);

  /// Extend the curve into the past so that it can be evaluated
  /// at these times. The times should be less than 
  virtual void extendBack(const std::vector<Time>& times,
                          const std::vector<Eigen::VectorXd>& values);


  /// \brief Fit the curve to these data points
  ///
  /// Underneath the curve should have some default policy for fitting
  virtual void fitCurve(const std::vector<int64_t>& times,
                        const std::vector<Eigen::VectorXd>& values);

  ///@}

  /// \name Methods to evaluate the curve.
  ///@{

  /// Evaluate the ambient space of the curve.
  virtual Eigen::VectorXd evaluateVector(Time time);
  
  /// Evaluate the curve derivatives.
  virtual Eigen::VectorXd evaluateDerivative(Time time, unsigned derivativeOrder);

  /// \brief Get an evaluator at this time
  EvaluatorTypePtr getTypedEvaluator(Time time);

  ///@}

 private:
  HermiteCoefficientManager manager_;
};

} // namespace curves


#endif /* CURVES_LINEAR_INTERPOLATION_VECTOR_SPACE_CURVE_HPP */
