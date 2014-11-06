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

namespace gtsam {
namespace traits {

template<>
struct equals<SE3> {
  typedef SE3 type;
  typedef bool result_type;
  bool operator()(const SE3& a, const SE3& b, double tol) {
    return a == b;
  }
};

template<>
struct print<SE3> {
  typedef SE3 type;
  typedef void result_type;
  void operator()(const SE3& obj, const std::string& str) {
    // todo replace print dummy
    std::cout << str << " " << std::endl;
  }
};
} // namespace traits
} // namespace gtsam

namespace gtsam {
typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
template<>
struct DefaultChart<SE3> {
  typedef SE3 type;
  typedef Eigen::Matrix<double, 6, 1> vector;
  static vector local(const SE3& origin, const SE3& other) {
    vector *diff;
    SE3CoefficientImplementation::localCoordinates(origin, other, diff);
    Eigen::Map<vector> map(diff->data());
    return vector(map);
  }
  static SE3 retract(const SE3& origin, const vector& d) {
    SE3 *retracted;
    SE3CoefficientImplementation::retract(origin, d, retracted);
    return *retracted;
  }
  static int getDimension(const SE3& origin) {
    return 6;
  }
};

} //namespace gtsam

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
  curve.print("slerp curve");

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
