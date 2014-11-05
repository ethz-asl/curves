#include <curves/SlerpSE3Curve.hpp>
#include <curves/SE3CoefficientImplementation.hpp>
#include <iostream>
#include "gtsam_unstable/nonlinear/Expression.h"

namespace curves {

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

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
  std::cout << "curve defined between times: " << manager_.getMinTime() << " and " << manager_.getMaxTime() <<std::endl;
  std::cout <<"=========================================" <<std::endl;
  for (size_t i = 0; i < manager_.size(); i++) {
    ss << "coefficient " << keys[i] << ": ";
    manager_.getCoefficientByKey(keys[i]).print(ss.str());
    std::cout << " | time: " << times[i];
    std::cout << std::endl;
    ss.str("");
  }
  std::cout <<"=========================================" <<std::endl;
}

void SlerpSE3Curve::getCoefficientsAt(const Time& time,
                                      Coefficient::Map* outCoefficients) const {
  CHECK_NOTNULL(outCoefficients);
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
  (*outCoefficients)[rval0->key] = rval0->coefficient;
  (*outCoefficients)[rval1->key] = rval1->coefficient;
}

void SlerpSE3Curve::getCoefficientsAt(const Time& time,
                                      KeyCoefficientTime** outCoefficient0,
                                      KeyCoefficientTime** outCoefficient1) const {
  CHECK_NOTNULL(&outCoefficient0);
  CHECK_NOTNULL(&outCoefficient1);
  if (time == this->getMaxTime()) {
    std::cout <<"max time reached" <<std::endl;
  }
  bool success = manager_.getCoefficientsAt(time, outCoefficient0, outCoefficient1);
  CHECK(success) << "Unable to get the coefficients at time " << time;
}

void SlerpSE3Curve::getCoefficientsInRange(Time startTime,
                                           Time endTime,
                                           Coefficient::Map* outCoefficients) const {
  manager_.getCoefficientsInRange(startTime, endTime, outCoefficients);
}

void SlerpSE3Curve::getCoefficients(Coefficient::Map* outCoefficients) const {
  manager_.getCoefficients(outCoefficients);
}

void SlerpSE3Curve::setCoefficient(Key key, const Coefficient& value) {
  manager_.setCoefficientByKey(key, value);
}

void SlerpSE3Curve::setCoefficients(const Coefficient::Map& coefficients) {
  manager_.setCoefficients(coefficients);
}

Time SlerpSE3Curve::getMaxTime() const {
  return manager_.getMaxTime();
}

Time SlerpSE3Curve::getMinTime() const {
  return manager_.getMinTime();
}

void SlerpSE3Curve::fitCurve(const std::vector<Time>& times,
                             const std::vector<ValueType>& values) {
  CHECK_EQ(times.size(), values.size());

  if(times.size() > 0) {
    Eigen::VectorXd val(7);
    manager_.clear();
    std::vector<Key> outKeys;
    outKeys.reserve(times.size());
    std::vector<Coefficient> coefficients;
    coefficients.reserve(times.size());
    for(size_t i = 0; i < values.size(); ++i) {
      CoefficientImplementation::Ptr impl(new SE3CoefficientImplementation);
      boost::dynamic_pointer_cast<SE3CoefficientImplementation>(impl)->makeValue(values[i],&val);
      Coefficient c1(impl,val);
      coefficients.push_back(c1);
    }
    manager_.insertCoefficients(times,coefficients,&outKeys);
  }
}

void SlerpSE3Curve::extend(const std::vector<Time>& times,
                           const std::vector<ValueType>& values) {

  //  CHECK_EQ(times.size(), values.size()) << "number of times and number of coefficients don't match";
  //  std::vector<Key> outKeys;
  //  std::vector<Coefficient> coefficients(values.size());
  //  for (size_t i = 0; i < values.size(); ++i) {
  //    coefficients[i] = Coefficient(values[i]);
  //  }
  //  manager_.insertCoefficients(times, coefficients, &outKeys);
}

typename SlerpSE3Curve::ValueType SlerpSE3Curve::evaluate(Time time) const {

  // \todo Jenkins implement this please
}

typename SlerpSE3Curve::DerivativeType SlerpSE3Curve::evaluateDerivative(Time time, unsigned derivativeOrder) const {

  // time is out of bound --> error
  CHECK_GE(time, this->getMinTime()) << "Time out of bounds";
  CHECK_LE(time, this->getMaxTime()) << "Time out of bounds";

  Eigen::VectorXd dCoeff;
  Time dt;
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);
  // first derivative
  if (derivativeOrder == 1) {
    dCoeff = rval1->coefficient.getValue() - rval0->coefficient.getValue();
    dt = rval1->time - rval0->time;
    return dCoeff/dt;
  } else { // order of derivative > 1 returns vector of zeros
    const int dimension = rval0->coefficient.dim();
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

// Evaluation function in functional form. To be passed to the expression
// old type Eigen::Matrix<double,dim,1> <double,dim,dim> <Eigen::MatrixXd*>&

SE3 evalFunc(SE3  v1, SE3  v2, double alpha,
             boost::optional<Eigen::Matrix<double,6,6>&>H1=boost::none,
             boost::optional<Eigen::Matrix<double,6,6>&>H2=boost::none,
             boost::optional<Eigen::Matrix<double,6,1>&>H3=boost::none) {

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

  // Rotational component computed without Baker-Campbell-Hausdorff
  //  a_AA_b.setAngle(-alpha_*a_AA_b.angle()/2);

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

  //TODO check matrix sizes should be chainRule.rows() x coefficient.ndim()
  if (H1) { *H1 += Eigen::Matrix<double,6,6>::Identity()*J_A_T; }
  if (H2) { *H2 += Eigen::Matrix<double,6,6>::Identity()*J_B_T; }
  if (H3) { *H3 = Eigen::Matrix<double,6,1>::Zero(); }

  // from previous evaluate function:

  SO3 a_R_b = w_R_a.inverted()*w_R_b;
  AngleAxis delta(a_R_b);
  delta.setAngle( delta.angle()*alpha);
  SE3 w_T_i(w_R_a*SO3(delta)  , (w_t_a*(1-alpha)+w_t_b*alpha).eval());

  //  return (w_T_i).getTransformationMatrix();
  return (w_T_i);
}

gtsam::Expression<typename SlerpSE3Curve::ValueType> SlerpSE3Curve::getEvalExpression(const Time& time) const {
  typedef typename SlerpSE3Curve::ValueType ValueType;
  using namespace gtsam;
  KeyCoefficientTime *rval0, *rval1;
  bool success = manager_.getCoefficientsAt(time, &rval0, &rval1);

  Expression<ValueType> leaf1(rval0->key);
  Expression<ValueType> leaf2(rval1->key);

  double alpha = double(time - rval0->time)/double(rval1->time - rval0->time);

  Expression<ValueType> rval(evalFunc, leaf1, leaf2, Expression<double>(alpha));

  return rval;
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


} // namespace curves
