/*
 * @file CompositionCurve.hpp
 * @date Feb 06, 2015
 * @author Abel Gawel, Renaud Dubé, Mike Bosse
 */

#ifndef CURVES_COMPOSITION_CURVE_HPP
#define CURVES_COMPOSITION_CURVE_HPP

namespace curves {

// SE3CompositionCurve is a curve composed of a base and a correction curve.
// The corrections can be sampled at a lower frequency than the base curve,
// Therefore reducing the optimization state space.
// The corrections are applied on the left side.

template <class C1, class C2>
class SE3CompositionCurve : public SE3Curve {
 private:
  C1 baseCurve_;
  C2 correctionCurve_;
 public:
  typedef SE3Curve::ValueType ValueType;
  typedef SE3Curve::DerivativeType DerivativeType;

  SE3CompositionCurve();
  ~SE3CompositionCurve();

  /// Print the value of the coefficient, for debugging and unit tests
    virtual void print(const std::string& str = "") const;

    /// The first valid time for the curve.
    virtual Time getMinTime() const;

    /// The one past the last valid time for the curve.
    virtual Time getMaxTime() const;

    bool isEmpty() const;

    // return number of coefficients curve is composed of
    int size() const;

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

    /// \brief Add coefficients to the correction curve at given times
    void addCorrectionCoefficients(const std::vector<Time>& times);

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
    virtual gtsam::Expression<ValueType> getValueExpression(const Time& time) const;

    virtual gtsam::Expression<DerivativeType> getDerivativeExpression(const Time& time, unsigned derivativeOrder) const;

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

    /// Initialize a GTSAM values structure with the desired keys
    virtual void initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) const;

    /// Initialize a GTSAM values structure for all keys
    virtual void initializeGTSAMValues(gtsam::Values* values) const;

    // updates the relevant curve coefficients from the GTSAM values structure
    virtual void updateFromGTSAMValues(const gtsam::Values& values);

};

} // namespace curves


#include "SE3CompositionCurve-inl.hpp"

#endif /* CURVES_COMPOSITION_CURVE_HPP */
