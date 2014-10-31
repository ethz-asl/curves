#include <gtest/gtest.h>
#include <curves/LinearInterpolationVectorSpaceCurve.hpp>
#include "gtsam_unstable/nonlinear/Expression.h"

#define DIM 3

using namespace curves;
using namespace gtsam;
using namespace std;

typedef typename LinearInterpolationVectorSpaceCurve<DIM>::ValueType ValueType;

namespace gtsam {
namespace traits {
// todo fix template specialization for eigen types. see :
// https://forum.kde.org/viewtopic.php?f=74&t=121280
template<>
bool equals(const Vector3& a,
            const Vector3& b,
            double tol) {
  return (a-b).cwiseAbs().maxCoeff() < tol;
}

template<>
void print(const Vector3& obj, const std::string& str) {
  std::cout << str << " " << obj << std::endl;
}

template<>
bool equals(const double& a,
            const double& b,
            double tol) {
  if (a > b)
    return a-b < tol;
  else
    return b-a < tol;
}

template<>
void print(const double& obj, const std::string& str) {
  std::cout << str << " " << obj << std::endl;
}
}
}

TEST(CurvesTestSuite, testGetExpression) {

  LinearInterpolationVectorSpaceCurve<DIM> curve;

  const double t[] = {0, 10};
  const double evalTime = 5;

  std::vector<Time> times(t,t+2);
  std::vector<ValueType> values;
  values.push_back(ValueType(0,0,0));
  values.push_back(ValueType(2,2,2));

  curve.fitCurve(times, values);

  KeyCoefficientTime *rval0, *rval1;
  curve.getCoefficientsAt(evalTime, &rval0, &rval1);

  Expression<ValueType> expression = curve.getEvalExpression(evalTime);

  std::set<Key> keys = expression.keys();

  // Check keys
  ASSERT_EQ(*(keys.begin()), rval0->key);
  ASSERT_EQ(*(++(keys.begin())), rval1->key);
  //  cout << "rval0->key: " << rval0->key << "rval1->key: " << rval1->key << endl;

  Values gtsamValues;
  gtsamValues.insert(rval0->key, ValueType(rval0->coefficient.getValue()));
  gtsamValues.insert(rval1->key, ValueType(rval1->coefficient.getValue()));
  JacobianMap actualMap;
  Eigen::MatrixXd H = Eigen::Matrix3d::Zero();
  actualMap.insert(make_pair(rval0->key, H.block(0, 0, 3, 3)));
  actualMap.insert(make_pair(rval1->key, H.block(0, 0, 3, 3)));
  ValueType result = expression.value(gtsamValues, actualMap);

  // Check evaluation result
  ASSERT_EQ(result, Eigen::Vector3d(1,1,1));

}



