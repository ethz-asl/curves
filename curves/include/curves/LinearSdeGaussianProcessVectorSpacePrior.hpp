#ifndef CURVES_LINEAR_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP
#define CURVES_LINEAR_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP

#include "GaussianProcessVectorSpacePrior.hpp"
#include "HermiteCoefficientManager.hpp"

namespace curves {

class LinearSdeGaussianProcessVectorSpacePrior : public GaussianProcessVectorSpacePrior {
 public:
  typedef GaussianProcessVectorSpacePrior Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::EvaluatorType EvaluatorType;
  typedef Parent::EvaluatorTypePtr EvaluatorTypePtr;
  typedef HermiteCoefficientManager CurveCoefficientManagerType;

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
  /// returns error if time is out of bounds
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
  virtual void fitCurve(const std::vector<Time>& times, const std::vector<ValueType>& values) {
    CHECK(false) << "The values of a Gaussian Process prior based on an SDE cannot be set, they are determined functionally.";
  }
  boost::unordered_map<Key, KeyCoefficientTime> getKeyCoefficientTime() const {
    CHECK(false) << "Not implemented for a Gaussian Process prior based on an SDE.";
  }

  ///@}

  const Eigen::MatrixXd& getPowerSpectralDensityMatrix() const {return stationaryPowerSpectralDensity_;}
  const Eigen::MatrixXd& getInversePowerSpectralDensityMatrix() const {return invStationaryPowerSpectralDensity_;}

 private:

  /// Stationary power spectral density matrix belonging to the zero-mean white noise process in the lienar SDE
  /// --- see Q_C in equation X, Anderson et al. (TBD)
  const Eigen::MatrixXd stationaryPowerSpectralDensity_;
  const Eigen::MatrixXd invStationaryPowerSpectralDensity_;

  /// Initial Conditions
  Time initialTime_;
  Eigen::MatrixXd initialInverseCovariance_;

  /// Exogenous Inputs
  std::map<Time, Eigen::VectorXd> timeToExogenousInputValue_;

  /// Evaluation of the prior mean function at `keytimes'
  std::map<Time, boost::shared_ptr<Eigen::VectorXd> > keytimeToMeanEvaluation_;

  /// Evaluations of other commonly required vectors/matrices based on the drift/diffusion of the prior
  boost::unordered_map<std::pair<Time,Time>, boost::shared_ptr<Eigen::VectorXd> > keytimesToLiftedExogenousInput_;
  boost::unordered_map<std::pair<Time,Time>, boost::shared_ptr<Eigen::MatrixXd> > keytimesToStateTransitionMatrix_;
  boost::unordered_map<std::pair<Time,Time>, boost::shared_ptr<Eigen::MatrixXd> > keytimesToInverseLiftedCovarianceMatrix_;

};

} // namespace curves

#endif /* CURVES_LTI_SDE_GAUSSIAN_PROCESS_VECTOR_SPACE_PRIOR_HPP */
