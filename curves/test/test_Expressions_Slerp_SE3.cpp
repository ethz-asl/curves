/*
 * @file test_Expressions_Slerp_SE3.cpp
 * @date Nov 05, 2014
 * @author Abel Gawel, Renaud Dube
 */

#include <gtest/gtest.h>
#include <curves/SlerpSE3Curve.hpp>
#include <curves/SE3CoefficientImplementation.hpp>
#include "gtsam_unstable/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "gtsam_unstable/nonlinear/ExpressionFactor.h"
#include <gtsam/base/numericalDerivative.h>

#include <Eigen/Core>
#include <boost/assign/list_of.hpp>

#define DIM 6

using namespace curves;
using namespace gtsam;
using namespace std;

typedef typename SlerpSE3Curve::ValueType ValueType;
typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;

// definition of traits for SE3 (to act as gtsam values)
namespace gtsam {
namespace traits {

template<>
struct equals<SE3> {
  typedef SE3 type;
  typedef bool result_type;
  bool operator()(const SE3& a, const SE3& b, double tol) {
    Eigen::Matrix<double,6,1> delta;
    return (a.log() -b.log()).cwiseAbs().maxCoeff() < tol;
  }
};

template<>
struct print<SE3> {
  typedef SE3 type;
  typedef void result_type;
  void operator()(const SE3& obj, const std::string& str) {
    // make nicer layout
    std:: cout << str << std::endl;
    std:: cout << "Position: " << obj.getPosition() << std::endl;
    std:: cout << "Rotation: "<< obj.getRotation() << std::endl;
  }
};

template <>
struct is_manifold<ChartValue<ValueType> > : public boost::true_type {};
} // namespace traits
} // namespace gtsam

// definition of charts for SE3, to be used in gtsam framework
namespace gtsam {
typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
template<>
struct DefaultChart<SE3> {
  typedef SE3 type;
  typedef Eigen::Matrix<double, 6, 1> vector;
  static vector local(const SE3& origin, const SE3& other) {
    vector diff;
    SE3CoefficientImplementation::localCoordinates(origin, other, &diff);
    return diff;
  }
  static SE3 retract(const SE3& origin, const vector& d) {
    SE3 retracted;
    SE3CoefficientImplementation::retract(origin, d, &retracted);
    return retracted;
  }
  static int getDimension(const SE3& origin) {
    return 6;
  }
};

} //namespace gtsam


// Wrapper to enable numerical differentiation on SlerpSE3Evaluator::evaluate
class ExpressionValueWrapper {
 private:
  Expression<ValueType> expression_;
  const Values& initialCoefficients_;

 public:
  ExpressionValueWrapper(const Expression<ValueType>& expression,
                         const Values& initialCoefficients) : expression_(expression),
                         initialCoefficients_(initialCoefficients){

  }
  ~ExpressionValueWrapper() {}

  ChartValue<ValueType> evaluate(const ChartValue<ValueType>& x, int arg) {
    std::set<Key>::const_iterator it = expression_.keys().begin();
    if (arg) {++it;}
    Key key = *it;

    Values values(initialCoefficients_);
    values.update(key,x.value());

    std::vector<size_t> dimensions;
    dimensions.push_back(DIM);
    dimensions.push_back(DIM);
    static const int Dim = traits::dimension<ValueType>::value;
    VerticalBlockMatrix Ab(dimensions, Dim);
    FastVector<Key> keys = boost::assign::list_of(*(expression_.keys().begin()))(*(++expression_.keys().begin()));
    JacobianMap actualMap(keys,Ab);

    ValueType result = expression_.value(values, actualMap);

    return convertToChartValue<ValueType>(result);
  }
};


// Expression to calculate relative measurement between 2 SE3 values
// \todo AG move this to trajectory maintainer at some point (this replaces relative pose factor functionality)
ValueType relativeMeasurementExpression(const ValueType& interp1,
                                        const ValueType& interp2,
                                        boost::optional<Eigen::Matrix<double,6,6>&>H1=boost::none,
                                        boost::optional<Eigen::Matrix<double,6,6>&>H2=boost::none) {
  if (H1) { *H1 = - Eigen::Matrix<double,6,6>::Identity(); }
  if (H2) { *H2 = Eigen::Matrix<double,6,6>::Identity(); }

  return ValueType(gtsam::DefaultChart<ValueType>::local(interp1, interp2));
}

// test for correct keys and evaluation function in Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3ExpressionKeysAndEvaluation) {
  SlerpSE3Curve curve;
  const double t[] = {0, 10};
  const double evalTime = 5;

  // setup two SE3 objects
  ValueType poseA(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  ValueType poseB(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(2,2,2));

  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(poseA);
  values.push_back(poseB);

  // interpolate curve
  curve.fitCurve(times, values);

  KeyCoefficientTime *rval0, *rval1;
  curve.getCoefficientsAt(evalTime, &rval0, &rval1);

  // get expression at evaluation time
  Expression<ValueType> expression = curve.getEvalExpression(evalTime);
  std::set<Key> keys = expression.keys();

  // Check keys
  ASSERT_EQ(*(keys.begin()), rval0->key);
  ASSERT_EQ(*(++(keys.begin())), rval1->key);

  // get coefficients from curve
  ValueType val0(SO3(rval0->coefficient.getValue()(3),rval0->coefficient.getValue().segment<3>(4)),rval0->coefficient.getValue().head<3>());
  ValueType val1(SO3(rval1->coefficient.getValue()(3),rval1->coefficient.getValue().segment<3>(4)),rval1->coefficient.getValue().head<3>());

  // fill retrieved coefficients in gtsam values container
  Values gtsamValues;
  gtsamValues.insert(rval0->key, val0);
  gtsamValues.insert(rval1->key, val1);

  // initialize objects
  std::vector<size_t> dimensions;
  dimensions.push_back(DIM);
  dimensions.push_back(DIM);
  static const int Dim = traits::dimension<ValueType>::value;
  VerticalBlockMatrix Ab(dimensions, Dim);
  FastVector<Key> key = boost::assign::list_of(rval0->key)(rval1->key);
  JacobianMap actualMap(key,Ab);

  // read out SE3 object from values container
  ValueType result = expression.value(gtsamValues, actualMap);
  Eigen::Vector3d resultPos = result.getPosition();
  Eigen::Vector4d resultRot = result.getRotation().vector();

  // assert return values are as expected
  ASSERT_EQ(resultPos, Eigen::Vector3d(1,1,1));
  ASSERT_NEAR(resultRot(0),0.9238795325112867,1e-6);
  ASSERT_NEAR(resultRot(1),0.0,1e-6);
  ASSERT_NEAR(resultRot(2),-0.3826834323650897,1e-6);
  ASSERT_NEAR(resultRot(3),0.0,1e-6);
}

// test for correct keys and evaluation function in Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3CurveEvaluate) {
  SlerpSE3Curve curve;
  const double t[] = {0, 10};

  // setup two SE3 objects
  ValueType poseA(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  ValueType poseB(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(2,2,2));

  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(poseA);
  values.push_back(poseB);

  // interpolate curve
  curve.fitCurve(times, values);

  // read out SE3 object from values container
  ValueType result = curve.evaluate(5);
  Eigen::Vector3d resultPos = result.getPosition();
  Eigen::Vector4d resultRot = result.getRotation().vector();

  // assert return values are as expected
  ASSERT_EQ(resultPos, Eigen::Vector3d(1,1,1));
  ASSERT_NEAR(resultRot(0),0.9238795325112867,1e-6);
  ASSERT_NEAR(resultRot(1),0.0,1e-6);
  ASSERT_NEAR(resultRot(2),-0.3826834323650897,1e-6);
  ASSERT_NEAR(resultRot(3),0.0,1e-6);

  result = curve.evaluate(0);
  resultPos = result.getPosition();
  resultRot = result.getRotation().vector();

  ASSERT_EQ(resultPos, Eigen::Vector3d(0,0,0));
  ASSERT_NEAR(resultRot(0),1,1e-6);
  ASSERT_NEAR(resultRot(1),0,1e-6);
  ASSERT_NEAR(resultRot(2),0,1e-6);
  ASSERT_NEAR(resultRot(3),0,1e-6);

  result = curve.evaluate(10);
  resultPos = result.getPosition();
  resultRot = result.getRotation().vector();

  ASSERT_EQ(resultPos, Eigen::Vector3d(2,2,2));
  ASSERT_NEAR(resultRot(0),0.7071067811865476,1e-6);
  ASSERT_NEAR(resultRot(1),0,1e-6);
  ASSERT_NEAR(resultRot(2),-0.7071067811865476,1e-6);
  ASSERT_NEAR(resultRot(3),0,1e-6);

}


// compare 2 types of expressions for Slerp SE3 interpolation:
// 1. Full Jacobian derivation & interpolation within 1 Expression (classic approach)
// 2. Assembly of expression by atomic Jacobians and operations (composition, inverse, exponential)
TEST(CurvesTestSuite, compareEvalExpressions1and2) {
  SlerpSE3Curve curve;
  const double t[] = {0, 10};
  const double evalTime = 5;

  // setup two SE3 objects
  // ValueType poseA(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  ValueType poseA(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(1,1,1));
  ValueType poseB(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(2,2,2));


  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(poseA);
  values.push_back(poseB);

  // interpolate curve
  curve.fitCurve(times, values);

  KeyCoefficientTime *rval0, *rval1;
  curve.getCoefficientsAt(evalTime, &rval0, &rval1);

  // get expression at evaluation time
  Expression<ValueType> expression1 = curve.getEvalExpression(evalTime);
  Expression<ValueType> expression2 = curve.getEvalExpression2(evalTime);

  std::set<Key> keys1 = expression1.keys();
  std::set<Key> keys2 = expression2.keys();

  // Check keys1
  ASSERT_EQ(*(keys1.begin()), rval0->key);
  ASSERT_EQ(*(++(keys1.begin())), rval1->key);
  ASSERT_EQ(*(keys2.begin()), rval0->key);
  ASSERT_EQ(*(++(keys2.begin())), rval1->key);

  // get coefficients from curve
  ValueType val0(SO3(rval0->coefficient.getValue()(3),rval0->coefficient.getValue().segment<3>(4)),rval0->coefficient.getValue().head<3>());
  ValueType val1(SO3(rval1->coefficient.getValue()(3),rval1->coefficient.getValue().segment<3>(4)),rval1->coefficient.getValue().head<3>());

  // fill retrieved coefficients in gtsam values container
  Values gtsamValues, gtsamValues2;
  gtsamValues.insert(rval0->key, val0);
  gtsamValues.insert(rval1->key, val1);

  // initialize objects
  std::vector<size_t> dimensions;
  dimensions.push_back(DIM);
  dimensions.push_back(DIM);
  static const int Dim = traits::dimension<ValueType>::value;
  VerticalBlockMatrix Ab(dimensions, Dim);
  Ab.matrix().setZero();
  FastVector<Key> key = boost::assign::list_of(rval0->key)(rval1->key);
  JacobianMap actualMap(key,Ab);
  ValueType result2 = expression2.value(gtsamValues, actualMap);

  // Test the jacobians
  // Call analytical calculation of Jacobians (using BAD)
  ExpressionValueWrapper expressionValueWrapper(expression2,gtsamValues);
  // Call numerical calculation of Jacobians
  gtsam::Matrix expectedH0 = gtsam::numericalDerivative11<ChartValue<ValueType>, ChartValue<ValueType> >
  (boost::bind(&ExpressionValueWrapper::evaluate, &expressionValueWrapper, _1, 0), convertToChartValue<ValueType>(val0), 1e-3);
  gtsam::Matrix expectedH1 = gtsam::numericalDerivative11<ChartValue<ValueType>, ChartValue<ValueType> >
  (boost::bind(&ExpressionValueWrapper::evaluate, &expressionValueWrapper, _1, 1), convertToChartValue<ValueType>(val1), 1e-3);

  // assert equality of Jacobians
  gtsam::Matrix analytical(actualMap(rval0->key));
  for (int i=0; i<analytical.size(); ++i){
    ASSERT_NEAR(expectedH0(i), analytical(i), 1e-3);
  }
}

// test basic gtsam interface of Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3ExpressionGTSAMoptimization) {

  SlerpSE3Curve curve;

  // Populate coefficients
  const double t[] = {0, 45, 90};
  std::vector<Time> times(t,t+3);
  std::vector<ValueType> coefficients;
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  coefficients.push_back(ValueType(SO3(0.9238795325112867,SO3::Vector3(0,-0.3826834323650897,0)),SE3::Position(4.5,4.5,4.5)));
  coefficients.push_back(ValueType(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(9,9,9)));

  // Populate measurements
  const double tmeas[] = {20, 40, 60, 80};
  std::vector<Time> measTimes(tmeas,tmeas+4);
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(SO3(0.984807753012208, SO3::Vector3(0,-0.17364817766693036,0)),SE3::Position(2,2,2)));
  measurements.push_back(ValueType(SO3(0.9396926207859083,SO3::Vector3(0,-0.3420201433256687,0)),SE3::Position(4,4,4)));
  measurements.push_back(ValueType(SO3(0.8660254037844386,SO3::Vector3(0,-0.5,0)),SE3::Position(6,6,6)));
  measurements.push_back(ValueType(SO3(0.7660444431189781,SO3::Vector3(0,-0.6427876096865393,0)),SE3::Position(8,8,8)));

  // Fit a curve
  curve.fitCurve(times, coefficients);

  // Populate GTSAM values
  KeyCoefficientTime *rval0, *rval1;
  Values initials, expected;
  for(int i=0; i< coefficients.size(); i++) {
    curve.getCoefficientsAt(times[i], &rval0, &rval1);
    Key key;
    if (i == coefficients.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    initials.insert(key,ValueType(SO3(1,SO3::Vector3(0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0)));
    expected.insert(key,coefficients[i]);
  }

  // Noise models
  noiseModel::Diagonal::shared_ptr measNoiseModel = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1));
  noiseModel::Diagonal::shared_ptr priorNoiseModel = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0));

  // Get expressions and build the graph
  NonlinearFactorGraph graph;
  for(int i=0; i < measurements.size(); i++) {

    Expression<ChartValue<ValueType> > predicted(convertToChartValue<ValueType>,
                                                 curve.getEvalExpression(measTimes[i]));

    ExpressionFactor<ChartValue<ValueType> > f(measNoiseModel,
                                               ChartValue<ValueType>(measurements[i]),
                                               predicted);
    graph.add(f);

    // Assert that error is null for expected values
    std::vector<gtsam::Matrix> H(2);
    Vector error = f.unwhitenedError(expected);
    CHECK((error).cwiseAbs().maxCoeff() < 1e-6);
  }
  // Add a prior of 0 on the first coefficient
  curve.getCoefficientsAt(0, &rval0, &rval1);
  Expression<ValueType> leaf1(rval0->key);
  Expression<ChartValue<ValueType> > prior(convertToChartValue<ValueType>, leaf1);
  graph.add(ExpressionFactor<ChartValue<ValueType> >(priorNoiseModel,ChartValue<ValueType>(ValueType(SO3(1.0,SO3::Vector3(0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0))),prior));

  // Optimize
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials).optimize();
  ASSERT_TRUE(expected.equals(result,1e-6));
}

// test for dynamic handling of duplicate keys in gtsam::ExpressionFactor
TEST(CurvesTestSuite, testSlerpSE3ExpressionDynamicKeys){
  SlerpSE3Curve curve;

  // Populate coefficients
  const double t[] = {0, 45, 90, 135};
  std::vector<Time> times(t,t+4);
  std::vector<ValueType> coefficients;
  coefficients.push_back(ValueType(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
  coefficients.push_back(ValueType(SO3(0.9238795325112867,SO3::Vector3(0,-0.3826834323650897,0)),SE3::Position(4.5,4.5,4.5)));
  coefficients.push_back(ValueType(SO3(0.7071067811865476,SO3::Vector3(0,-0.7071067811865476,0)),SE3::Position(9,9,9)));
  coefficients.push_back(ValueType(SO3(0.38268343236508984,SO3::Vector3(0,-0.9238795325112866,0)),SE3::Position(13.5,13.5,13.5)));

  // Populate measurements
  const double tmeas[] = {10, 30, 20};
  std::vector<Time> measTimes(tmeas,tmeas+3);
  const double durations[] = {30, 30, 100};
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(SO3(0.9659258262890683,SO3::Vector3(0,-0.25881904510252074,0)),SE3::Position(3.633974596215561,3,2.633974596215561)));
  measurements.push_back(ValueType(SO3(0.9659258262890683,SO3::Vector3(0,-0.25881904510252074,0)),SE3::Position(4.901923788646684,3,1.901923788646684)));
  measurements.push_back(ValueType(SO3(0.6427876096865394,SO3::Vector3(0,-0.766044443118978,0)),SE3::Position(14.316911861358276,10,10.377680849309444)));

  // Fit a curve
  curve.fitCurve(times, coefficients);
  noiseModel::Diagonal::shared_ptr measNoiseModel = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1));

  // Populate GTSAM values
  KeyCoefficientTime *rval0, *rval1;
  Values initials, expected;
  for(int i=0; i< coefficients.size(); i++) {
    curve.getCoefficientsAt(times[i], &rval0, &rval1);
    Key key;
    if (i == coefficients.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    expected.insert(key,coefficients[i]);
  }

  // create ExpressionFactor which represents relative pose relation (former Relative pose factor)
  for (int i=0; i<3; ++i) {
    Expression<ChartValue<ValueType> > relative (relativeMeasurementExpression, curve.getEvalExpression(tmeas[i]), curve.getEvalExpression(tmeas[i]+durations[i]));
    ExpressionFactor<ChartValue<ValueType> > factor(measNoiseModel,ChartValue<ValueType>(measurements[i]), relative);
    const FastVector<Key> keys = factor.keys();

    // check that the ExpressionFactor dynamically handles duplicate keys
    ASSERT_EQ(factor.keys().size(),i+2);
    Values::iterator itGTSAM = expected.begin();
    for(FastVector<Key>::const_iterator it = factor.keys().begin(); it!= factor.keys().end(); ++it) {
      ASSERT_EQ(*it,(*itGTSAM).key);
      ++itGTSAM;
    }
    Vector error = factor.unwhitenedError(expected);
    CHECK((error).cwiseAbs().maxCoeff() < 1e-6);
  }
}

TEST(CurvesTestSuite, testSlerpSE3ExpressionJacobians){
  // todo AG-RD
  ASSERT_TRUE(true);
}
