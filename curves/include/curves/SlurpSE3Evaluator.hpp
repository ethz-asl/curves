#ifndef CURVES_SLURP_SE3_EVALUATOR_HPP
#define CURVES_SLURP_SE3_EVALUATOR_HPP

#include "Evaluator.hpp"
#include "SE3Config.hpp"
#include "SlurpSE3Curve.hpp"

namespace curves {

class Coefficients;

class SlurpSE3Evaluator : public Evaluator<SE3Config> {

 public:

  typedef Evaluator<SE3Config> Parent;
  typedef Parent::DerivativeType DerivativeType;
  typedef Parent::ValueType ValueType;

  SlurpSE3Evaluator(const SlurpSE3Curve& curve, const Time& time);
  virtual ~SlurpSE3Evaluator();

  virtual void getKeys(std::vector<Key> *outKeys) const;

  virtual void appendKeys(std::vector<Key> *outKeys) const;

  virtual std::vector<Key>::const_iterator keyBegin() const;

  virtual std::vector<Key>::const_iterator keyEnd() const;

  /// Evaluate the ambient space of the curve (functional form) by specifying new coefficients.
  virtual ValueType evaluate(const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the curve derivatives (functional form).
  virtual ValueType evaluateDerivative(unsigned derivativeOrder,
                                             const std::vector<Coefficient>& coefficients) const;

  /// Evaluate the ambient space of the curve (functional form).
  virtual ValueType evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                                    std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the curve derivatives (functional form).
  virtual ValueType evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                        const std::vector<Coefficient>& coefficients,
                                                        std::vector<Eigen::MatrixXd>* outJacobian) const;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluate(const Coefficients& coefficients) const;

  /// Get the curve Jacobians.
  /// This is the main interface for GTSAM
  virtual void getJacobians(unsigned derivativeOrder,
                            const Coefficients& coefficients,
                            const Eigen::MatrixXd& chainRule,
                            const std::vector<Eigen::MatrixXd*>& jacobians) const;

  /// Evaluate the ambient space of the curve
  virtual ValueType evaluateDerivative(unsigned derivativeOrder,
                                       const Coefficients& coefficients) const;
 private:

  std::vector<Key> keys_;
  double alpha_;
  size_t dimension_;

};

} // namespace curves


#endif /* CURVES_SLURP_SE3_EVALUATOR_HPP */
