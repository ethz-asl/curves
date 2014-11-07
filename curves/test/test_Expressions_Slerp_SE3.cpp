#include <gtest/gtest.h>
#include <curves/SlerpSE3Curve.hpp>
#include <curves/SE3CoefficientImplementation.hpp>
#include "gtsam_unstable/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "gtsam_unstable/nonlinear/ExpressionFactor.h"
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

// Expression to calculate relative measurement between 2 SE3 values
// \todo AG move this to trajectory maintainer at some point (this replaces relative pose factor functionality)
ValueType relativeMeasurementExpression(const ValueType& interp1,
                                        const ValueType& interp2,
                                        boost::optional<Eigen::Matrix<double,6,6>&>H1=boost::none,
                                        boost::optional<Eigen::Matrix<double,6,6>&>H2=boost::none) {
  if (H1) { *H1 = Eigen::Matrix<double,6,6>::Identity(); }
  if (H2) { *H2 = - Eigen::Matrix<double,6,6>::Identity(); }

  return ValueType(gtsam::DefaultChart<ValueType>::local(interp1, interp2));
}

// test for correct keys and evaluation function in Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3ExpressionKeysAndEvaluation) {

  SlerpSE3Curve curve;
  const double t[] = {0, 10};
  const double evalTime = 5;

  ValueType poseA(SO3(SO3::Vector4(1,0,0,0)),SE3::Position(0,0,0));
  ValueType poseB(SO3(SO3::Vector4(0.7071067811865476,0,-0.7071067811865476,0)),SE3::Position(2,2,2));

  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(poseA);
  values.push_back(poseB);

  curve.fitCurve(times, values);

  KeyCoefficientTime *rval0, *rval1;
  curve.getCoefficientsAt(evalTime, &rval0, &rval1);

  Expression<ValueType> expression = curve.getEvalExpression(evalTime);

  std::set<Key> keys = expression.keys();

  // Check keys
  ASSERT_EQ(*(keys.begin()), rval0->key);
  ASSERT_EQ(*(++(keys.begin())), rval1->key);

  Values gtsamValues;

  gtsamValues.insert(rval0->key, ValueType(SO3(SO3::Vector4(rval0->coefficient.getValue().segment<4>(3))),rval0->coefficient.getValue().head<3>()));
  gtsamValues.insert(rval1->key, ValueType(SO3(SO3::Vector4(rval1->coefficient.getValue().segment<4>(3))),rval1->coefficient.getValue().head<3>()));

  Eigen::MatrixXd H = Eigen::Matrix3d::Zero();
  std::vector<size_t> dimensions;
  dimensions.push_back(DIM);
  dimensions.push_back(DIM);
  static const int Dim = traits::dimension<ValueType>::value;
  VerticalBlockMatrix Ab(dimensions, Dim);
  FastVector<Key> key = boost::assign::list_of(rval0->key)(rval1->key);
  JacobianMap actualMap(key,Ab);
  ValueType result = expression.value(gtsamValues, actualMap);
  Eigen::Vector3d resultPos = result.getPosition();
  Eigen::Vector4d resultRot = result.getRotation().vector();

  ASSERT_EQ(resultPos, Eigen::Vector3d(1,1,1));
  ASSERT_NEAR(resultRot(0),0.9238795325112867,1e-6);
  ASSERT_NEAR(resultRot(1),0.0,1e-6);
  ASSERT_NEAR(resultRot(2),-0.3826834323650897,1e-6);
  ASSERT_NEAR(resultRot(3),0.0,1e-6);
}

// test basic gtsam interface of Slerp SE3 curves
TEST(CurvesTestSuite, testSlerpSE3ExpressionGTSAMoptimization) {

  SlerpSE3Curve curve;

  // Populate coefficients
  const double t[] = {0, 45, 90};
  std::vector<Time> times(t,t+3);
  std::vector<ValueType> coefficients;

  coefficients.push_back(ValueType(SO3(SO3::Vector4(1,0,0,0)),SE3::Position(0,0,0)));
  coefficients.push_back(ValueType(SO3(SO3::Vector4(0.9238795325112867,0,-0.3826834323650897,0)),SE3::Position(4.5,4.5,4.5)));
  coefficients.push_back(ValueType(SO3(SO3::Vector4(0.7071067811865476,0,-0.7071067811865476,0)),SE3::Position(9,9,9)));

  // Populate measurements
  const double tmeas[] = {20, 40, 60, 80};
  std::vector<Time> measTimes(tmeas,tmeas+4);
  std::vector<ValueType> measurements;

  measurements.push_back(ValueType(SO3(SO3::Vector4(0.984807753012208,0,-0.17364817766693036,0)),SE3::Position(2,2,2)));
  measurements.push_back(ValueType(SO3(SO3::Vector4(0.9396926207859083,0,-0.3420201433256687,0)),SE3::Position(4,4,4)));
  measurements.push_back(ValueType(SO3(SO3::Vector4(0.8660254037844386,0,-0.5,0)),SE3::Position(6,6,6)));
  measurements.push_back(ValueType(SO3(SO3::Vector4(0.7660444431189781,0,-0.6427876096865393,0)),SE3::Position(8,8,8)));

  // Fit a curve
  curve.fitCurve(times, coefficients);

  // Populate GTSAM values
  KeyCoefficientTime *rval0, *rval1;
  Values initials, expected;
  for(int i=0; i< coefficients.size(); i++) {
    curve.getCoefficientsAt(times[i], &rval0, &rval1);
    Key key;
    // todo use another method for getting the keys
    // here a if is necessary since t=maxtime the coef is in rval1
    if (i == coefficients.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    initials.insert(key,ValueType(SO3(SO3::Vector4(1.0,0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0)));
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
  graph.add(ExpressionFactor<ChartValue<ValueType> >(priorNoiseModel,ChartValue<ValueType>(ValueType(SO3(SO3::Vector4(1.0,0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0))),prior));

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

  // create a set of coefficients
  coefficients.push_back(ValueType(SO3(SO3::Vector4(1,0,0,0)),SE3::Position(0,0,0)));
  coefficients.push_back(ValueType(SO3(SO3::Vector4(0.9238795325112867,0,-0.3826834323650897,0)),SE3::Position(4.5,4.5,4.5)));
  coefficients.push_back(ValueType(SO3(SO3::Vector4(0.7071067811865476,0,-0.7071067811865476,0)),SE3::Position(9,9,9)));
  coefficients.push_back(ValueType(SO3(SO3::Vector4(0.38268343236508984,0,-0.9238795325112866,0)),SE3::Position(13.5,13.5,13.5)));

  // Populate measurements
  const double tmeas[] = {10, 30, 20};
  std::vector<Time> measTimes(tmeas,tmeas+3);
  const double durations[] = {30, 30, 100};
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(SO3(SO3::Vector4(0.9659258262890683,0,-0.25881904510252074,0)),SE3::Position(3,3,3)));
  measurements.push_back(ValueType(SO3(SO3::Vector4(0.9659258262890683,0,-0.25881904510252074,0)),SE3::Position(3,3,3)));
  measurements.push_back(ValueType(SO3(SO3::Vector4(0.6427876096865394,0,-0.766044443118978,0)),SE3::Position(10,10,10)));

  // Fit a curve
  curve.fitCurve(times, coefficients);
  noiseModel::Diagonal::shared_ptr measNoiseModel = noiseModel::Diagonal::Sigmas((gtsam::Vector(6) << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1));

  // create ExpressionFactor which represents relative pose relation (former Relative pose factor)
  for (int i=0; i<3; ++i) {
    Expression<ChartValue<ValueType> > relative (relativeMeasurementExpression, curve.getEvalExpression(tmeas[i]), curve.getEvalExpression(tmeas[i]+durations[i]));
    ExpressionFactor<ChartValue<ValueType> > factor(measNoiseModel,ChartValue<ValueType>(measurements[i]), relative);
    const FastVector<Key> keys = factor.keys();
    // check that the ExpressionFactor dynamically handles duplicate keys
    ASSERT_EQ(factor.keys().size(),i+2);
  }
}

TEST(CurvesTestSuite, testSlerpSE3ExpressionJacobians){
  // todo AG-RD
  ASSERT_TRUE(true);
}
