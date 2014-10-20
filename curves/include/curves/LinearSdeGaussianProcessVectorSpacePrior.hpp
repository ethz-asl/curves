#ifndef CURVES_LINEAR_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_LINEAR_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "GaussianProcessVectorSpacePrior.hpp"
#include "HermiteCoefficientManager.hpp"

namespace curves {

/// \class LinearSdeGaussianProcessVectorSpacePrior
///
/// Interface for vector-space, Gaussian-process-trajectory priors based on
/// linear stochastic different equations.
class LinearSdeGaussianProcessVectorSpacePrior : public GaussianProcessVectorSpacePrior {
 public:
  /// \brief Parent class
  typedef GaussianProcessVectorSpacePrior Parent;

  /// \brief The value type of the curve.
  typedef Parent::ValueType ValueType;

  /// \brief The derivative type of the curve.
  typedef Parent::DerivativeType DerivativeType;

  /// \brief The evaluator type of the curve.
  typedef Parent::EvaluatorType EvaluatorType;

  /// \brief The evaluator type pointer.
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;

  /// \brief The coefficient manager type the GP curve should use.
  typedef HermiteCoefficientManager CurveCoefficientManagerType;

  /// \struct LinearSdeCoefficient
  ///
  /// Stores useful variable evaluations at and between keytimes.
  struct LinearSdeCoefficient {
    /// Timestamp of the prior evaluation
    Time time;

    /// Evaluation of the prior's mean
    Eigen::VectorXd mean;

    /// Previous evaluation's timestamp
    Time prevTime;

    /// Integrated exogenous input between the current
    /// and previous evaluation timestamps
    Eigen::VectorXd liftedExogenousInput;

    /// State transition matrix evaluation between the
    /// current and previous evaluation timestamps
    Eigen::MatrixXd stateTransitionMatrix;

    /// Integrated white noise and diffusion between the
    /// current and previous evaluation timestamps
    Eigen::MatrixXd inverseLiftedCovarianceMatrix;
  };

  /// \brief Initialize with the dimension of the vector space
  LinearSdeGaussianProcessVectorSpacePrior(size_t dimension, Eigen::MatrixXd stationaryPowerSpectralDensity);
  virtual ~LinearSdeGaussianProcessVectorSpacePrior();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// The first valid time for the curve.
  virtual Time getMinTime() const;

  /// The one past the last valid time for the curve.
  virtual Time getMaxTime() const;

  /// Initialize the prior
  virtual void initialize(Time initialTime, Eigen::VectorXd initialMean, Eigen::MatrixXd initialInverseCovariance);
  virtual bool isInitialized() const;

  /// Append an exogenous input.
  /// Exogenous inputs create a step function that is assumed to hold constant from last known value.
  virtual void addExogenousInput(Time time, const ValueType& value);

  /// Set the discrete time exogenous inputs.
  /// Exogenous inputs create a step function that is assumed to hold constant from last known value.
  virtual void setExogenousInputs(const std::vector<Time>& times,
                                  const std::vector<ValueType>& values);

  /// Get the number of keytimes.
  virtual unsigned getNumKeyTimes() const;


  /// Evaluate the prior mean.
  virtual Eigen::VectorXd evaluate(Time time) const;

  /// Evaluate the derivative of the prior mean.
  virtual Eigen::VectorXd evaluateDerivative(Time time, unsigned derivativeOrder) const;

  /// Evaluate the prior at the query time and the key times (associated
  /// with local support), and evaluate the interpolation matrix
  /// associated with the key times (non-zero blocks of K(t)K^{-1}).
  virtual Eigen::VectorXd evaluateAndInterpMatrices(Time time, const std::vector<Time>& keyTimes,
                                                    std::vector<Eigen::VectorXd>* outEvalAtKeyTimes,
                                                    std::vector<Eigen::MatrixXd>* outInterpMatrices) const;

  /// Evaluate the derivative of the prior at the query time and the
  /// key times (associated with local support), and evaluate the
  /// interpolation matrix associated with the derivative order and
  /// key times (non-zero blocks of K(t)K^{-1}).
  virtual Eigen::VectorXd evaluateDerivativeAndInterpMatrices(Time time, unsigned derivativeOrder, const std::vector<Time>& keyTimes,
                                                              std::vector<Eigen::VectorXd>* outEvalAtKeyTimes,
                                                              std::vector<Eigen::MatrixXd>* outInterpMatrices) const;

  /// \brief Get an evaluator at this time
  EvaluatorTypePtr getEvaluator(const Time& time) const;

  /// \brief Set time range of prior
  virtual void setTimeRange(Time minTime, Time maxTime);

  /// To be implemented by SDE-form-specific class
  /// --- see v_i in equation X, Anderson et al. (TBD)
  /// \todo needs better name
  virtual Eigen::VectorXd calculateLiftedExogenousInput(Time time1, Time time2) const = 0;

  /// To be implemented by SDE-form-specific class
  /// --- see Phi(t1,t2) in equation X, Anderson et al. (TBD)
  virtual Eigen::MatrixXd calculateStateTransitionMatrix(Time time1, Time time2) const = 0;

  /// To be implemented by SDE-form-specific class
  /// --- see Q_i in equation X, Anderson et al. (TBD)
  /// \todo needs better name
  virtual Eigen::MatrixXd calculateLiftedCovarianceMatrix(Time time1, Time time2) const = 0;
  virtual Eigen::MatrixXd calculateInverseLiftedCovarianceMatrix(Time time1, Time time2) const = 0;

  /// \name Methods `removed' from curve functionality
  /// \todo better way to do this?
  ///@{

  /// \brief Unused for GP priors based on SDEs
  virtual void getCoefficientsAt(const Time& time, Coefficient::Map* outCoefficients) const {
    CHECK(false) << "Not implemented for a Gaussian Process prior based on an SDE.";
  }
  virtual void getCoefficientsInRange(Time startTime, Time endTime, Coefficient::Map* outCoefficients) const {
    CHECK(false) << "Not implemented for a Gaussian Process prior based on an SDE.";
  }
  virtual void getCoefficients(Coefficient::Map* outCoefficients) const {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }
  virtual void setCoefficient(Key key, const Coefficient& value) {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }
  virtual void setCoefficients(const Coefficient::Map& coefficients){
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }
  virtual void extend(const std::vector<Time>& times, const std::vector<ValueType>& values) {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }
  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values, std::vector<Key>* outKeys = NULL) {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }
  ///@}

  /// Get the power spectral density matrix
  const Eigen::MatrixXd& getPowerSpectralDensityMatrix() const {return stationaryPowerSpectralDensity_;}

  /// Get the inverse of the power spectral density matrix
  const Eigen::MatrixXd& getInversePowerSpectralDensityMatrix() const {return invStationaryPowerSpectralDensity_;}

 private:

  /// Stationary power spectral density matrix belonging to the zero-mean white noise process in the linear SDE
  /// --- see Q_C in equation X, Anderson et al. (TBD)
  const Eigen::MatrixXd stationaryPowerSpectralDensity_;
  const Eigen::MatrixXd invStationaryPowerSpectralDensity_;

  /// Initial Conditions
  bool initialized_;
  Time initialTime_;
  Eigen::VectorXd initialMean_;
  Eigen::MatrixXd initialCovariance_;

  /// Exogenous Inputs
  std::map<Time, Eigen::VectorXd> timeToExogenousInputValue_;

  /// Coefficient maps (ordered and unordered for use cases)
  std::map<Time, boost::shared_ptr<LinearSdeCoefficient> > keytimeToMeanSorted_;
  boost::unordered_map<Time, boost::shared_ptr<LinearSdeCoefficient> > keytimeToMean_;

  /// Update the mean function due to a change in the exogenous
  /// input at some time. Note in the case of multiple changes,
  /// one call to the earliest time should suffice.
  void updateFromExogenousInputChange(Time time);

  /// Add an exogenous input change to the prior.
  void addExogenousInput(Time time, const ValueType& value, bool updateMean);

  /// Add a keytime to the prior
  virtual void addKeyTime(const Time& time);

  /// Add a keytimes to the prior
  virtual void addKeyTimes(const std::vector<Time>& times);

  /// Clear the keytimes in the prior
  virtual void clearKeyTimes();
};

} // namespace curves

#endif /* CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
