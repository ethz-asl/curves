/*
 * @file SlerpSE3Curve.cpp
 * @date Oct 10, 2014
 * @author Renaud Dube, Abel Gawel
 */

#include <curves/SlerpSE3Curve.hpp>
#include <iostream>

namespace curves {

SlerpSE3Curve::SlerpSE3Curve() : SE3Curve() {}

SlerpSE3Curve::~SlerpSE3Curve() {}

void SlerpSE3Curve::print(const std::string& str) const {
  std::cout << "=========================================" << std::endl;
  std::cout << "=========== Slerp SE3 CURVE =============" << std::endl;
  std::cout << str << std::endl;
  std::cout << "num of coefficients: " << manager_.size() << std::endl;
  std::cout << "dimension: " << 6 << std::endl;
  std::stringstream ss;
  std::vector<Key> keys;
  std::vector<Time> times;
  manager_.getTimes(&times);
  manager_.getKeys(&keys);
  std::cout << "curve defined between times: " << manager_.getMinTime() <<
      " and " << manager_.getMaxTime() <<std::endl;
  std::cout <<"=========================================" <<std::endl;
  for (size_t i = 0; i < manager_.size(); i++) {
    ss << "coefficient " << keys[i] << ": ";
    gtsam::traits<Coefficient>::Print(manager_.getCoefficientByKey(keys[i]),ss.str());
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
}

Time SlerpSE3Curve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time SlerpSE3Curve::getMinTime() const {
  return manager_.getMinTime();
}

bool SlerpSE3Curve::isEmpty() const {
  std::vector<Time> outTimes;
  manager_.getTimes(&outTimes);
  return outTimes.empty();
}

int SlerpSE3Curve::size() const {
  return manager_.size();
}

void SlerpSE3Curve::fitCurve(const std::vector<Time>& times,
                             const std::vector<ValueType>& values,
                             std::vector<Key>* outKeys) {
  CHECK_EQ(times.size(), values.size());
  if(times.size() > 0) {
    manager_.insertCoefficients(times,values, outKeys);
  }
}

void SlerpSE3Curve::extend(const std::vector<Time>& times,
                           const std::vector<ValueType>& values) {

  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  std::vector<Key> outKeys;

  manager_.insertCoefficients(times, values, &outKeys);
}

typename SlerpSE3Curve::DerivativeType
SlerpSE3Curve::evaluateDerivative(Time time,
                                  unsigned derivativeOrder) const {

  // time is out of bound --> error
  CHECK_GE(time, this->getMinTime()) << "Time out of bounds";
  CHECK_LE(time, this->getMaxTime()) << "Time out of bounds";

  Eigen::VectorXd dCoeff;
  Time dt;
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  // first derivative
  if (derivativeOrder == 1) {
    //todo Verify this
    dCoeff = gtsam::traits<Coefficient>::Local(rval1->coefficient,rval0->coefficient);
    dt = rval1->time - rval0->time;
    return dCoeff/dt;
  } else { // order of derivative > 1 returns vector of zeros
    const int dimension = gtsam::traits<Coefficient>::dimension;
    return Eigen::VectorXd::Zero(dimension,1);
  }
}

Eigen::Matrix3d crossOperator(Eigen::Vector3d vector){
  Eigen::Matrix3d rval;
  rval << 0.0, -vector(2), vector(1),
      vector(2), 0.0, -vector(0),
      -vector(1), vector(0), 0.0;
  return rval;
}

/// \brief Evaluation function in functional form. To be passed to the expression (full implementation)
SE3 slerpInterpolation(SE3  v1, SE3  v2, double alpha,
                       gtsam::OptionalJacobian<6,6> H1,
                       gtsam::OptionalJacobian<6,6> H2) {

  typedef Eigen::Matrix3d Matrix3d;
  SO3 w_R_a(v1.getRotation());
  SO3 w_R_b(v2.getRotation());
  SE3::Position w_t_a(v1.getPosition());
  SE3::Position w_t_b(v2.getPosition());

  w_R_a.getRotationMatrix();
  w_R_b.getRotationMatrix();

  Matrix3d J_tA_tT, J_rA_tT, J_tA_rT, J_rA_rT;
  Matrix3d J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
  Eigen::Matrix<double, 6, 6> J_A_T, J_B_T;

  AngleAxis a_AA_b(w_R_a.inverted() * w_R_b);

  // Rotational component computed with Baker-Campbell-Hausdorff
  Matrix3d Re = alpha *(Matrix3d::Identity() - ((1-alpha)/2)*a_AA_b.angle()*crossOperator(a_AA_b.axis()));

  J_rA_tT = -((1-alpha) * crossOperator(w_t_a) + alpha*crossOperator(w_t_b))*Re + alpha*crossOperator(w_t_b);
  J_tA_rT = Matrix3d::Zero();
  J_rA_rT = Matrix3d::Identity() - Re;
  J_tA_tT = Matrix3d::Identity()*(1.0 - alpha);

  J_A_T << J_tA_tT, J_rA_tT, J_tA_rT, J_rA_rT;

  J_tB_rT = Matrix3d::Zero();
  J_rB_rT = Re;
  J_tB_tT = Matrix3d::Identity()*alpha;
  J_rB_tT = (crossOperator(w_t_a)*(1-alpha)+crossOperator(w_t_b)*alpha)*Re - alpha*crossOperator(w_t_b);

  J_B_T << J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;

  /// \todo check matrix sizes should be chainRule.rows() x coefficient.ndim()
  if (H1) { *H1 += J_A_T; }
  if (H2) { *H2 += J_B_T; }

  /// Slerp interpolation function
  SO3 a_R_b = w_R_a.inverted()*w_R_b;
  AngleAxis delta(a_R_b);
  delta.setAngle( delta.angle()*alpha);
  SE3 w_T_i(w_R_a*SO3(delta)  , (w_t_a*(1-alpha)+w_t_b*alpha).eval());

  return (w_T_i);
}

/// \brief \f[T^{\alpha}\f]
SE3 transformationPower(SE3  T, double alpha,
                        gtsam::OptionalJacobian<6,6> H) {
  typedef Eigen::Matrix3d Matrix3d;

  SO3 R(T.getRotation());
  SE3::Position t(T.getPosition());

  AngleAxis angleAxis(R);
  angleAxis.setUnique();
  angleAxis.setAngle( angleAxis.angle()*alpha);
  angleAxis.setUnique();

  if(H) {
    Matrix3d J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
    Eigen::Matrix<double, 6, 6> J_B_T;

    AngleAxis AA(R);

    Matrix3d Re = alpha *(Matrix3d::Identity() - ((1-alpha)/2)*AA.angle()*crossOperator(AA.axis()));

    J_tB_rT = Matrix3d::Zero();
    J_rB_rT = Re;
    J_tB_tT = Matrix3d::Identity()*alpha;
    J_rB_tT = (crossOperator(t)*alpha)*Re - alpha*crossOperator(t);

    J_B_T << J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
    (*H) = J_B_T;
  }

  return SE3(SO3(angleAxis),(t*alpha).eval());
}

/// \brief \f[A*B\f]
SE3 composeTransformations(SE3 A, SE3 B,
                           gtsam::OptionalJacobian<6,6> H1,
                           gtsam::OptionalJacobian<6,6> H2) {
  // todo compute jacobians
  if (H1) {
    (*H1) = Eigen::Matrix<double,6,6>::Identity();
  }

  if (H2) {
    Eigen::Matrix3d J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
    Eigen::Matrix<double, 6, 6> J_B_T;
    SE3::Position ta(A.getPosition());
    J_tB_rT = Eigen::Matrix3d::Zero();
    J_rB_rT = A.getRotationMatrix();
    J_tB_tT = A.getRotationMatrix();
    J_rB_tT = crossOperator(ta)*A.getRotationMatrix();

    J_B_T << J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
    (*H2) = J_B_T;
  }

  return A*B;
}

/// \brief \f[T^{-1}\f]
SE3 inverseTransformation(SE3 T, gtsam::OptionalJacobian<6,6> H) {
  // todo compute jacobian
  if (H) {
    Eigen::Matrix3d J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
    Eigen::Matrix<double, 6, 6> J_B_T;
    SE3::Position ta(T.getPosition());
    J_tB_rT = Eigen::Matrix3d::Zero();
    J_rB_rT = - T.getRotationMatrix().inverse();
    J_tB_tT = - T.getRotationMatrix().inverse();
    J_rB_tT = T.getRotationMatrix().inverse() * crossOperator(ta);

    J_B_T << J_tB_tT, J_rB_tT, J_tB_rT, J_rB_rT;
    (*H) = J_B_T;
  }

  return T.inverted();
}

Eigen::Vector3d transformPoint(SE3 A, Eigen::Vector3d p, gtsam::OptionalJacobian<3,6> H) {

  Eigen::Vector3d Tp = A * p;
  if (H) {
    Eigen::Matrix3d J_tA, J_rA;
    Eigen::Matrix<double, 3, 6> J_A;
    J_tA = Eigen::Matrix3d::Identity();
    //todo: is J_tp_tT maybe Identity?
    J_rA = -crossOperator(Tp);

    J_A << J_tA, J_rA;
    (*H) = J_A;
  }
  return Tp;
}

/// \brief forms slerp interpolation into a binary expression with 2 leafs and binds alpha into it
///        \f[ T = A(A^{-1}B)^{\alpha} \f]
gtsam::Expression<typename SlerpSE3Curve::ValueType>
SlerpSE3Curve::getEvalExpression(const Time& time) const {
  typedef typename SlerpSE3Curve::ValueType ValueType;
  using namespace gtsam;
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);

  Expression<ValueType> leaf1(rval0->key);
  Expression<ValueType> leaf2(rval1->key);

  double alpha = double(time - rval0->time)/double(rval1->time - rval0->time);

  Expression<ValueType> rval(boost::bind(&slerpInterpolation,_1,_2,alpha,_3,_4),
                             leaf1, leaf2);

  return rval;
}
/// \brief forms slerp interpolation into a binary expression with 2 leafs and binds alpha into it,
///        uses break down of expression into its operations
///        \f[ T = A(A^{-1}B)^{\alpha} \f]
gtsam::Expression<typename SlerpSE3Curve::ValueType>
SlerpSE3Curve::getEvalExpression2(const Time& time) const {
  typedef typename SlerpSE3Curve::ValueType ValueType;
  using namespace gtsam;
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);

  double alpha = double(time - rval0->time)/double(rval1->time - rval0->time);

  Expression<ValueType> leaf1(rval0->key);
  Expression<ValueType> leaf2(rval1->key);
  //  std::cout << "alpha" << alpha << std::endl;
  if (alpha == 0) {
    return leaf1;
  } else if (alpha == 1) {
    return leaf2;
  } else {
    Expression<ValueType> inverted(inverseTransformation, leaf1);
    Expression<ValueType> composed(composeTransformations, inverted, leaf2);
    Expression<ValueType> powered(boost::bind(&transformationPower,_1,alpha,_2), composed);
    return Expression<ValueType>(composeTransformations, leaf1, powered);
  }
}


SE3 SlerpSE3Curve::evaluate(Time time) const {
  KeyCoefficientTime *a, *b;
  bool success = manager_.getCoefficientsAt(time, &a, &b);
  double alpha = double(time - a->time)/double(b->time - a->time);

  SE3 delta = composeTransformations(inverseTransformation(a->coefficient),b->coefficient);

  return composeTransformations(pose0, transformationPower(delta,alpha));
}


void SlerpSE3Curve::setTimeRange(Time minTime, Time maxTime) {
  // \todo Abel and Renaud
  CHECK(false) << "Not implemented";
}

/// \brief Evaluate the angular velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the velocity of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearVelocityB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistA(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateTwistB(Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the angular derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateAngularDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame b as seen from Frame a, expressed in Frame a.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief Evaluate the derivative of Frame a as seen from Frame b, expressed in Frame b.
Eigen::Vector3d SlerpSE3Curve::evaluateLinearDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular derivative of Frame b as seen from Frame a,
/// expressed in Frame a. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeA(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}
/// \brief evaluate the velocity/angular velocity of Frame a as seen from Frame b,
/// expressed in Frame b. The return value has the linear velocity (0,1,2),
/// and the angular velocity (3,4,5).
Vector6d SlerpSE3Curve::evaluateDerivativeB(unsigned derivativeOrder, Time time) {
  CHECK(false) << "Not implemented";
}

void SlerpSE3Curve::initializeGTSAMValues(gtsam::FastVector<gtsam::Key> keys, gtsam::Values* values) {
  for (unsigned int i = 0; i < keys.size(); ++i) {
    values->insert(keys[i],manager_.getCoefficientByKey(keys[i]));
  }
}

void SlerpSE3Curve::updateFromGTSAMValues(const gtsam::Values& values) {
  gtsam::Values::const_iterator iter;
  for (iter = values.begin(); iter != values.end(); ++iter) {
    manager_.setCoefficientByKey(iter->key,iter->value.cast<Coefficient>());
  }
}

} // namespace curves
