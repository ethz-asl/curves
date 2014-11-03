#include <gtest/gtest.h>
#include <curves/LinearInterpolationVectorSpaceCurve.hpp>
#include "gtsam_unstable/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "gtsam_unstable/nonlinear/ExpressionFactor.h"

#include <boost/assign/list_of.hpp>

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

TEST(CurvesTestSuite, testExpressionKeysAndEvaluation) {

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

  Values gtsamValues;
  gtsamValues.insert(rval0->key, ValueType(rval0->coefficient.getValue()));
  gtsamValues.insert(rval1->key, ValueType(rval1->coefficient.getValue()));

  Eigen::MatrixXd H = Eigen::Matrix3d::Zero();
  std::vector<size_t> dimensions;
  dimensions.push_back(DIM);
  dimensions.push_back(DIM);
  static const int Dim = traits::dimension<ValueType>::value;
  VerticalBlockMatrix Ab(dimensions, Dim);
  FastVector<Key> key = boost::assign::list_of(rval0->key)(rval1->key);
  JacobianMap actualMap(key,Ab);
  ValueType result = expression.value(gtsamValues, actualMap);

  // Check evaluation result
  ASSERT_EQ(result, Eigen::Vector3d(1,1,1));
}

TEST(CurvesTestSuite, testExpressionGTSAMoptimization) {

  LinearInterpolationVectorSpaceCurve<DIM> curve;

  // Populate coefficients
  const double t[] = {0, 10, 20};
  const double evalTime = 5;
  std::vector<Time> times(t,t+3);
  std::vector<ValueType> coefficients;
  coefficients.push_back(ValueType(0,0,0));
  coefficients.push_back(ValueType(40,20,-30));
  coefficients.push_back(ValueType(20,-20,-80));

  // Populate measurements
  const double tmeas[] = {4, 8, 12, 16};
  std::vector<Time> measTimes(tmeas,tmeas+4);
  std::vector<ValueType> measurements;
  measurements.push_back(ValueType(16,8,-12));
  measurements.push_back(ValueType(32,16,-24));
  measurements.push_back(ValueType(36,12,-40));
  measurements.push_back(ValueType(28,-4,-60));

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
    initials.insert(key,ValueType(0,0,0));
    expected.insert(key,coefficients[i]);
  }

  // Noise models
  SharedNoiseModel measNoiseModel = noiseModel::Diagonal::Sigmas(Vector3(1, 1, 1));
  SharedNoiseModel priorNoiseModel = noiseModel::Diagonal::Sigmas(Vector3(0, 0, 0));

  // Get expressions and build the graph
  NonlinearFactorGraph graph;
  for(int i=0; i < measurements.size(); i++) {
    Expression<ChartValue<ValueType> > predicted(convertToChartValue<ValueType>,
                                                 curve.getEvalExpression(measTimes[i]));

    ExpressionFactor<ChartValue<ValueType> > f(measNoiseModel,
                                                ChartValue<Vector3>(measurements[i]),
                                                predicted);
    graph.add(ExpressionFactor<ChartValue<ValueType> >(measNoiseModel,
                                                       ChartValue<Vector3>(measurements[i]),
                                                       predicted));
    // Assert that error is null for expected values
    std::vector<Matrix> H(2);
    Vector error = f.unwhitenedError(expected, H);
    ASSERT_EQ(error, ValueType(0,0,0));
  }

  // Add a prior of 0,0,0 on the first coefficient
  curve.getCoefficientsAt(0, &rval0, &rval1);
  Expression<ValueType> leaf1(rval0->key);
  Expression<ChartValue<ValueType> > prior(convertToChartValue<ValueType>, leaf1);
  graph.add(ExpressionFactor<ChartValue<ValueType> >(priorNoiseModel,ChartValue<ValueType>(ValueType(0,0,0)),prior));

  // Optimize
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials).optimize();

  ASSERT_TRUE(expected.equals(result,1e-8));
}

TEST(CurvesTestSuite, testExpressionJacobians){
  // todo AG-RD
  ASSERT_TRUE(true);
}




