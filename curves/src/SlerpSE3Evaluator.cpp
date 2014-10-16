#include <curves/SlerpSE3Evaluator.hpp>
#include <curves/Coefficients.hpp>
#include <iostream>
#include "kindr/minimal/quat-transformation.h"

#include <unsupported/Eigen/MatrixFunctions>

using namespace std;

Eigen::Matrix3d crossOperator(Eigen::Vector3d vector){
  Eigen::Matrix3d rval;
  rval << 0.0, -vector(2), vector(1),
      vector(2), 0.0, -vector(0),
      -vector(1), vector(0), 0.0;
  return rval;
}

namespace curves {

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

SlerpSE3Evaluator::SlerpSE3Evaluator(const SlerpSE3Curve& curve, const Time& time){
  KeyCoefficientTime *coeff0, *coeff1;
  curve.getCoefficientsAt(time, &coeff0, &coeff1);
  keys_.push_back(coeff0->key);
  keys_.push_back(coeff1->key);

  // Compute alpha_
  Time dt = coeff1->time - coeff0->time;
  if(dt == 0) {
    alpha_ = 1.0;
  } else {
    Time t = time - coeff0->time;
    alpha_ = double(t)/double(dt);
  }
}

SlerpSE3Evaluator::~SlerpSE3Evaluator() {}

void SlerpSE3Evaluator::getKeys(std::vector<Key> *outKeys) const {
  CHECK_NOTNULL(outKeys);
  *outKeys = keys_;
}

void SlerpSE3Evaluator::appendKeys(std::vector<Key> *outKeys) const {
  CHECK_NOTNULL(outKeys);
  outKeys->insert(outKeys->end(), keys_.begin(), keys_.end());
}

std::vector<Key>::const_iterator SlerpSE3Evaluator::keyBegin() const {
  return keys_.begin();
}

std::vector<Key>::const_iterator SlerpSE3Evaluator::keyEnd() const {
  return keys_.end();
}


SlerpSE3Evaluator::ValueType SlerpSE3Evaluator::evaluate(
    const std::vector<Coefficient>& coefficients) const {
//  const Eigen::VectorXd& coeffA = coefficients[0].getValue();
//  const Eigen::VectorXd& coeffB = coefficients[1].getValue();
//
//  SE3 w_T_a(SO3(SO3::Vector4(coeffA.segment<4>(3))),coeffA.head<3>());
//  SE3 w_T_b(SO3(SO3::Vector4(coeffB.segment<4>(3))),coeffB.head<3>());
//
//  SE3 a_T_b = w_T_a.inverted()*w_T_b;
//  SE3::Vector6 delta = a_T_b.log()*alpha_;
//
//  SE3 a_T_i = SE3(delta);
//
//  return (w_T_a*a_T_i).getTransformationMatrix();
  CHECK(false) << "Implement this if necessary.";
}

SlerpSE3Evaluator::ValueType SlerpSE3Evaluator::evaluate(
    const Coefficients& coefficients) const {
  const Eigen::VectorXd& coeffA = coefficients.get(keys_[0]).getValue();
  const Eigen::VectorXd& coeffB = coefficients.get(keys_[1]).getValue();

  SO3 w_R_a(SO3::Vector4(coeffA.segment<4>(3)));
  SO3 w_R_b(SO3::Vector4(coeffB.segment<4>(3)));

  SE3::Position w_t_a(coeffA.head<3>());
  SE3::Position w_t_b(coeffB.head<3>());

  SO3 a_R_b = w_R_a.inverted()*w_R_b;
  AngleAxis delta(a_R_b);
  delta.setAngle( delta.angle()*alpha_ );
  SE3 w_T_i(w_R_a*SO3(delta)  , (w_t_a*(1-alpha_)+w_t_b*alpha_).eval());

  return (w_T_i).getTransformationMatrix();
}

void SlerpSE3Evaluator::getJacobians(unsigned derivativeOrder,
                                                           const Coefficients& coefficients,
                                                           const Eigen::MatrixXd& chainRule,
                                                           const std::vector<Eigen::MatrixXd*>& jacobians) const {

  typedef Eigen::Matrix3d Matrix3d;

  // TODO(Abel and Renaud) implement velocity
  CHECK_EQ(derivativeOrder, 0);
  CHECK_NOTNULL(jacobians[0]);
  CHECK_NOTNULL(jacobians[1]);
  CHECK_EQ(6,chainRule.rows());
  CHECK_EQ(6,chainRule.cols());

  const Eigen::VectorXd& coeffA = coefficients.get(keys_[0]).getValue();
  const Eigen::VectorXd& coeffB = coefficients.get(keys_[1]).getValue();
  SO3 w_R_a(SO3::Vector4(coeffA.segment<4>(3)));
  SO3 w_R_b(SO3::Vector4(coeffB.segment<4>(3)));
  SE3::Position w_t_a(coeffA.head<3>());
  SE3::Position w_t_b(coeffB.head<3>());

  Matrix3d Ra = w_R_a.getRotationMatrix();
  Matrix3d Rb = w_R_b.getRotationMatrix();

  Matrix3d J_tA_tT, J_rA_tT, J_tA_rT, J_rA_rT;
  Matrix3d J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
  Eigen::Matrix<double, 6, 6> J_A_T, J_B_T;



  Matrix3d Re = alpha_ * Ra * (Ra.transpose() * Rb).pow(-1/2) * Ra.transpose();

  J_rA_tT = -((1-alpha_) * crossOperator(w_t_a) + alpha_*crossOperator(w_t_b))*Re + alpha_*crossOperator(w_t_b);
  J_tA_rT = Matrix3d::Zero();
  J_rA_rT = Matrix3d::Identity() - Re;
  J_tA_tT = Matrix3d::Identity()*(1.0 - alpha_);

  J_A_T << J_tA_tT, J_rA_tT, J_tA_rT, J_rA_rT;

  J_tB_rT = Matrix3d::Zero();
  J_rB_rT = Re;
  J_tB_tT = Matrix3d::Identity()*alpha_;
  J_rB_tT = (crossOperator(w_t_a)*(1-alpha_)+crossOperator(w_t_b)*alpha_)*Re - alpha_*crossOperator(w_t_b);

  J_B_T << J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;

  std::cout << "chainrule: " << std::endl << chainRule << std::endl;

  //TODO check matrix sizes should be chainRule.rows() x coefficient.ndim()
  *(jacobians[0]) +=  chainRule*J_A_T;
  *(jacobians[1]) +=  chainRule*J_B_T;
}

/// Evaluate the ambient space of the curve
SlerpSE3Evaluator::ValueType SlerpSE3Evaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                            const Coefficients& coefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives (functional form).
SlerpSE3Evaluator::ValueType SlerpSE3Evaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                            const std::vector<Coefficient>& coefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the ambient space of the curve (functional form).
SlerpSE3Evaluator::ValueType SlerpSE3Evaluator::evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                                                             std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives (functional form).
SlerpSE3Evaluator::ValueType SlerpSE3Evaluator::evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                                                       const std::vector<Coefficient>& coefficients,
                                                                                       std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

} // namespace curves
