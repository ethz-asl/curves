/*
 * @file CompositionCurve.hpp
 * @date Feb 06, 2015
 * @author Renaud Dubé, Abel Gawel, Mike Bosse
 */

#include "curves/SE3Curve.hpp"

#pragma once

namespace curves {

// SE3CompositionCurve is a curve composed of a base and a correction curve.
// The corrections can be sampled at a lower frequency than the base curve,
// therefore reducing the optimization state space. The corrections are applied
// on the left side.

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

    /// \brief Print the value of the base and corrections curves coefficients
    virtual void print(const std::string& str = "") const;

    /// \brief Save curves data as .csv file
    void saveCurves(const std::string& filename) const;

    /// \brief Returns the first valid time for the curve.
    virtual Time getMinTime() const;

    /// \brief Returns the last valid time for the curve.
    virtual Time getMaxTime() const;

    /// \brief Checks if the curve is empty.
    bool isEmpty() const;

    /// \brief Returns the number of coefficients in the correction curve
    int size() const;

    /// \brief Returns the number of coefficients in the base curve
    int baseSize() const;

    /// \brief Returns the number of coefficients in the correction curve
    int correctionSize() const;

    /// \brief Extend the curve so that it can be evaluated at these times by
    ///        using a default correction sampling policy.
    virtual void extend(const std::vector<Time>& times,
                        const std::vector<ValueType>& values,
                        std::vector<Key>* outKeys = NULL);

    /// \brief Set the minimum sampling period for the correction curve.
    ///        Overloads the function defined in SE3Curve base class.
    void setMinSamplingPeriod(const Time minSamplingPeriod);

    /// \brief Set the sampling ratio for the correction curve.
    ///   eg. 4 will add a coefficient every 4 extend
    void setSamplingRatio(const int ratio);

    /// \brief Fold in the correction curve into the base curve and reinitialize
    ///        correction curve coefficients to identity transformations.
    void foldInCorrections();

    /// \brief Fit a new curve to these data points.
    ///
    /// The existing curve will be cleared.
    /// Underneath the curve should have some default policy for fitting.
    virtual void fitCurve(const std::vector<Time>& times,
                          const std::vector<ValueType>& values,
                          std::vector<Key>* outKeys = NULL);

    /// \brief Add coefficients to the correction curve at given times.
    void setCorrectionTimes(const std::vector<Time>& times);

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

    /// \brief Clear the base and correction curves.
    virtual void clear();

    /// \brief Remove a correction coefficient at the specified time.
    void removeCorrectionCoefficientAtTime(Time time);

    /// \brief Set the correction coefficient value at the specified time.
    void setCorrectionCoefficientAtTime(Time time, ValueType value);

    /// \brief Perform a rigid transformation on the left side of the curve
    void transformCurve(const ValueType T);

    /// \brief Reset the correction curve to identity values with knots at desired times
    void resetCorrectionCurve(const std::vector<Time>& times);

    /// \brief Set the base curve to given values with knots at desired times
    /// Resets the curve beforehand.
    void setBaseCurve(const std::vector<Time>& times, const std::vector<ValueType>& values);

    /// \brief Add / replace the given coefficients without resetting the curve.
    void setBaseCurvePart(const std::vector<Time>& times, const std::vector<ValueType>& values);

    /// \brief Modifies values of the base coefficient in batch, starting at times[0] and assuming that
    /// a coefficient exists at all the specified times.
    void modifyBaseCoefficientsValuesInBatch(const std::vector<Time>& times, const std::vector<ValueType>& values);

    /// \brief Save the base curve times and composed curve values
    void saveCurveTimesAndValues(const std::string& filename) const;

    void saveCurveAtTimes(const std::string& filename, std::vector<Time> times) const;

    void saveCorrectionCurveAtTimes(const std::string& filename, std::vector<Time> times) const;

    void saveCorrectionCurveTimesAndValues(const std::string& filename) const;

    void getBaseCurveTimes(std::vector<Time>* outTimes) const;

    void getBaseCurveTimesInWindow(std::vector<Time>* outTimes, Time begTime, Time endTime) const;

    void getCurveTimes(std::vector<Time>* outTimes) const;
};

} // namespace curves

#include "SE3CompositionCurve-inl.hpp"
