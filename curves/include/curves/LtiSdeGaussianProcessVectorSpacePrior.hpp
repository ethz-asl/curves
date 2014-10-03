#ifndef CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "GaussianProcessVectorSpacePrior.hpp"
#include "HermiteCoefficientManager.hpp"

namespace curves {

class LtiSdeGaussianProcessVectorSpacePrior : public GaussianProcessVectorSpacePrior {
 public:
  typedef GaussianProcessVectorSpacePrior Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::EvaluatorType EvaluatorType;
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;
  typedef HermiteCoefficientManager CurveCoefficientManagerType;

  /// \brief Initialize with the dimension of the vector space
  //LtiSdeGaussianProcessVectorSpacePrior(size_t dimension);
  LtiSdeGaussianProcessVectorSpacePrior(Eigen::MatrixXd ltiSdeMatrixF, Eigen::MatrixXd ltiSdeMatrixL, Eigen::MatrixXd ltiSdeMatrixQc);
  virtual ~LtiSdeGaussianProcessVectorSpacePrior();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// \brief Get the coefficients that are active at a certain time.
  virtual void getCoefficientsAt(const Time& time,
                                 Coefficient::Map* outCoefficients) const;

  /// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
  virtual void getCoefficientsInRange(Time startTime, 
                                      Time endTime, 
                                      Coefficient::Map* outCoefficients) const;

  /// \brief Get all of the curve's coefficients.
  virtual void getCoefficients(Coefficient::Map* outCoefficients) const;

  /// The first valid time for the curve.
  virtual Time getMinTime() const;

  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const;

  /// Append an exogenous input.
  /// Exogenous inputs create a step function that is assumed to hold constant from last known value.
  virtual void addExogenousInput(Time time, const ValueType& value);

  /// Set the discrete time exogenous inputs.
  /// Exogenous inputs create a step function that is assumed to hold constant from last known value.
  virtual void setExogenousInputs(const std::vector<Time>& times,
                                  const std::vector<ValueType>& values);

  virtual void addKeyTime(const Time& time);

  virtual void addKeyTimes(const std::vector<Time>& times);

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

  /// Evaluate the prior at the query time, the key times (associated with local support), and evaluate the interpolation matrix associated with the key times (and belonging to K(t)K^{-1}).
  virtual Eigen::VectorXd evaluateAndInterpMatrices(Time time, const std::vector<Time>& keyTimes,
                                                    const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                    const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const;

  virtual Eigen::VectorXd evaluateDerivativeAndInterpMatrices(Time time, unsigned derivativeOrder, const std::vector<Time>& keyTimes,
                                                              const std::vector<Eigen::VectorXd*>& outEvalAtKeyTimes,
                                                              const std::vector<Eigen::MatrixXd*>& outInterpMatrices) const;

  /// \brief Get an evaluator at this time
  EvaluatorTypePtr getEvaluator(const Time& time) const;

  virtual void setTimeRange(Time minTime, Time maxTime);

  /// Returns the Key-Coefficient-Time-relationship
  boost::unordered_map<Key, KeyCoefficientTime> getKeyCoefficientTime() const;

  /// \name Methods `removed' from curve functionality
  /// \todo better way to do this?
  ///@{

  /// \brief Unused for GP priors, a prior cannot be fit to values
  virtual void setCoefficient(Key key, const Coefficient& value) {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }

  /// \brief Unused for GP priors, a prior cannot be fit to values
  virtual void setCoefficients(const Coefficient::Map& coefficients){
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }

  /// \brief Unused for GP priors, a prior cannot be fit to values
  virtual void extend(const std::vector<Time>& times, const std::vector<ValueType>& values) {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }

  /// \brief Unused for GP priors, a prior cannot be fit to values
  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values) {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }

  ///@}

  // set analytical state transition function
  // --- should check it is a valid state transition matrix
  // --- should check if it is valid for the provided F

  // set analytical Q function
  // --- check?

  // set analytical Qinv function
  // --- check?

 private:
  HermiteCoefficientManager manager_;

  const Eigen::MatrixXd ltiSdeMatrixF_;
  const Eigen::MatrixXd ltiSdeMatrixL_;
  const Eigen::MatrixXd ltiSdeMatrixQc_;

  // IC val, time and cov
  // map<times, inputvals>

  // map<keytimes, shared<vals>]>

  // analytical state transition * = NULL
  // map<keytimes, shared<state transitions>]>

  // analytical Q * = NULL
  // map<keytime_pair, shared<Qis>]>

  /// Time to coefficient mapping
  //std::map<Time, Eigen::VectorXd> timeToValue_;
  //boost::unordered_map<Time, Eigen::MatrixXd> timeToStateTransitionMatrix_;
  //boost::unordered_map<std::pair<Time,Time>, Eigen::MatrixXd> timePairToQiMatrix_;

  // evalStateTransition(t1, t2)
  // evalQ(t1, t2)
  // evalQinv(t1, t2)

  //bool hasKeytimeAtTime(Time time, std::map<Time, Eigen::VectorXd>::iterator *it);


};

} // namespace curves

#endif /* CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
