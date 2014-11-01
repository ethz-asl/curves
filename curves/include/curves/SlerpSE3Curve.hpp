#ifndef CURVES_SLERP_SE3_CURVE_HPP
#define CURVES_SLERP_SE3_CURVE_HPP

#include "SE3Curve.hpp"
#include "HermiteCoefficientManager.hpp"

class SlerpSE3Evaluator; // Forward declaration

namespace curves {

class SlerpSE3Curve : public SE3Curve {
 public:
  typedef SE3Curve::ValueType ValueType;
  typedef SE3Curve::DerivativeType DerivativeType;
  typedef SE3Curve::EvaluatorType EvaluatorType;
  typedef SE3Curve::EvaluatorTypePtr EvaluatorTypePtr;

  SlerpSE3Curve();
  virtual ~SlerpSE3Curve();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// \brief Get the coefficients that are active at a certain time.
  virtual void getCoefficientsAt(const Time& time,
                                 Coefficient::Map* outCoefficients) const;

  /// \brief Get the KeyCoefficientTimes that are active at a certain time.
  void getCoefficientsAt(const Time& time,
                         KeyCoefficientTime** outCoefficient0,
                         KeyCoefficientTime** outCoefficient1) const;

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
                        const std::vector<ValueType>& values);

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

  /// \brief Get an evaluator at this time
  EvaluatorTypePtr getEvaluator(const Time& time) const;

  virtual gtsam::Expression<ValueType> getEvalExpression(const Time& time) const;

  virtual void setTimeRange(Time minTime, Time maxTime);

  /// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularVelocityA(Time time);

  /// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularVelocityB(Time time);

  /// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearVelocityA(Time time);

  /// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearVelocityB(Time time);

  /// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistA(Time time);

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateTwistB(Time time);

  /// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateAngularDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateAngularDerivativeB(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
  virtual Eigen::Vector3d evaluateLinearDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
  virtual Eigen::Vector3d evaluateLinearDerivativeB(unsigned derivativeOrder, Time time);

  /// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeA(unsigned derivativeOrder, Time time);

  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
  ///        and the angular velocity (3,4,5).
  virtual Vector6d evaluateDerivativeB(unsigned derivativeOrder, Time time);

 private:
  HermiteCoefficientManager manager_;
};

} // namespace curves

#include "SlerpSE3Evaluator.hpp"

#endif /* CURVES_SLERP_SE3_CURVE_HPP */
