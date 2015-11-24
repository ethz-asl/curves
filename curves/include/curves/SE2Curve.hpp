/*
 * @file SE2Curve.hpp
 * @date Nov 24, 2015
 * @author Renaud Dubé
 */

#ifndef SE2_CURVE_H_
#define SE2_CURVE_H_

#include "SE2Config.hpp"
#include "Curve.hpp"

namespace curves {

// Curves over SE2 inherit the interface from Curve and CurveBase and define specific
// methods to support physical interpretations of temporal derivatives.
//
// For the purposes of these uses, the curve is defined between Frame a and Frame b
// such that evaluate() returns \f$ \mathbf{T}_{a,b} \f$, the transformation that takes
// points from Frame b to Frame a.
//
class SE2Curve : public Curve<SE2Config> {
 public:
  SE2Curve();
  virtual ~SE2Curve();

  typedef Curve<SE2Config> Parent;
  typedef Parent::ValueType ValueType;
  typedef Parent::DerivativeType DerivativeType;

  //todo revisit these functions if needed
//  /// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
//  virtual Eigen::Vector3d evaluateAngularVelocityA(Time time) = 0;
//
//  /// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
//  virtual Eigen::Vector3d evaluateAngularVelocityB(Time time) = 0;
//
//  /// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
//  virtual Eigen::Vector3d evaluateLinearVelocityA(Time time) = 0;
//
//  /// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
//  virtual Eigen::Vector3d evaluateLinearVelocityB(Time time) = 0;
//
//  /// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
//  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
//  ///        and the angular velocity (3,4,5).
//  virtual Vector6d evaluateTwistA(Time time) = 0;
//
//  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
//  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
//  ///        and the angular velocity (3,4,5).
//  virtual Vector6d evaluateTwistB(Time time) = 0;
//
//  /// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
//  virtual Eigen::Vector3d evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) = 0;
//
//  /// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
//  virtual Eigen::Vector3d evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) = 0;
//
//  /// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
//  virtual Eigen::Vector3d evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) = 0;
//
//  /// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
//  virtual Eigen::Vector3d evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) = 0;
//
//  /// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
//  ///        expressed in Frame a. The return value has the linear velocity (0,1,2),
//  ///        and the angular velocity (3,4,5).
//  virtual Vector6d evaluateDerivativeA(unsigned derivativeOrder, Time time) = 0;
//
//  /// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
//  ///        expressed in Frame b. The return value has the linear velocity (0,1,2),
//  ///        and the angular velocity (3,4,5).
//  virtual Vector6d evaluateDerivativeB(unsigned derivativeOrder, Time time) = 0;

  /// \brief Get the dimension of this curve
  //virtual size_t dim() const;
  ///@}
 private:
};

}  // namespace curves

#endif // SE2_CURVE_H_
