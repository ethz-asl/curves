#ifndef CT_LINEAR_INTERPOLATION_CURVE_HPP
#define CT_LINEAR_INTERPOLATION_CURVE_HPP

namespace curves {

class LinearInterpolationCurve : public VectorSpaceCurve {
 public:
  LinearInterpolationCurve();
  virtual ~LinearInterpolationCurve();

    /// \brief Get the coefficients that are active at a certain time.
  virtual void getCoefficientsAt(int64_t time, 
                                 CoefficientMap& outCoefficients) const;

  /// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
  virtual void getCoefficientsInRange(int64_t startTime, 
                                      int64_t endTime, 
                                      CoefficientMap& outCoefficients) const;

  /// \brief Get all of the curve's coefficients.
  virtual void getCoefficients(CoefficientMap& outCoefficients) const;
  
  /// \brief Set a coefficient.
  virtual void setCoefficient(size_t key, const Coefficient& value);

  /// \brief Set coefficients.
  virtual void setCoefficients(CoefficientMap& coefficients);
  ///@}


  /// \name Methods to deal with time.
  ///@{

  /// The first valid time for the curve.
  virtual int64_t getBackTime() const;
  
  /// The one past the last valid time for the curve.
  virtual int64_t getFrontTime() const;

  ///@}

  /// \name Methods to grow or shrink the curve.

  /// Extend the curve into the future so that it can be evaluated
  /// at this time.
  virtual void extendFront(int64_t time);

  /// Extend the curve into the past so that it can be evaluated
  /// at this time.
  virtual void extendBack(int64_t time);

  /// Retract the curve in the front. The caller guarantees that
  /// they will never try to evaluate this time or higher again.
  virtual void retractFront(int64_t time);

  /// Retract the curve in the back. The caller guarantees that
  /// they will never try to evaluate this time or lower again.
  virtual void retractBack(int64_t time);

  /// Extend the curve into the future so that it can be evaluated
  /// at the following times. The arguments are smoothing points.
  virtual void extendFront(const std::vector<int64_t>& times,
                           const std::vector<Eigen::VectorXd>& values);

  /// Extend the curve into the past so that it can be evaluated
  /// at these times. The times should be less than 
  virtual void extendBack(const std::vector<int64_t>& times,
                          const std::vector<Eigen::VectorXd>& values);


  ///@}

  /// \name Methods to evaluate the curve.
  ///@{

  /// Evaluate the ambient space of the curve.
  virtual Eigen::VectorXd evaluateVector(int64_t time);
  
  /// Evaluate the curve derivatives.
  virtual Eigen::VectorXd evaluateDerivative(int64_t time, unsigned derivativeOrder);
};

} // namespace curves


#endif /* CT_LINEAR_INTERPOLATION_CURVE_HPP */
