/*
 * @file test_SE3_Sinus_Circle.cpp
 * @date Nov 05, 2014
 * @author Abel Gawel, Renaud Dube
 *
 *  Note : this test should eventually be transfered to trajectory maintainer
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

#include <iostream>
#include <fstream>
#include <stdlib.h>

#define DIM 6

using namespace curves;
using namespace std;
using namespace gtsam;

typedef SlerpSE3Curve::ValueType ValueType;
typedef SlerpSE3Curve::EvaluatorTypePtr EvaluatorTypePtr;

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

// definition of traits for SE3 (to act as gtsam values)
namespace gtsam {
namespace traits {

template<>
struct equals<SE3> {
  typedef SE3 type;
  typedef bool result_type;
  bool operator()(const SE3& a, const SE3& b, double tol) {
    Eigen::Matrix<double,6,1> delta;

    cout << "error " << (a.log() -b.log()).cwiseAbs().maxCoeff() << endl;
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

//// Expression to calculate relative measurement between 2 SE3 values
//// \todo AG move this to trajectory maintainer at some point (this replaces relative pose factor functionality)
//ValueType relativeMeasurementExpression(const ValueType& interp1,
//                                        const ValueType& interp2,
//                                        boost::optional<Eigen::Matrix<double,6,6>&>H1=boost::none,
//                                        boost::optional<Eigen::Matrix<double,6,6>&>H2=boost::none) {
//  if (H1) { *H1 = Eigen::Matrix<double,6,6>::Identity(); }
//  if (H2) { *H2 = - Eigen::Matrix<double,6,6>::Identity(); }
//
//  return ValueType(gtsam::DefaultChart<ValueType>::local(interp1, interp2));
//}


// test a greater example of using more factors, creating gtsam factor graph and optimizing it
TEST(CurvesTestSuite, test_MITb) {

  typedef SlerpSE3Curve::ValueType ValueType;
  typedef SlerpSE3Curve::EvaluatorTypePtr EvaluatorTypePtr;

  // read the relative measurements from file
  curves::Time timestamp;
  std::string value;
  std::ifstream relativeData ("MITb_edges_measurements.csv");
  CHECK(relativeData.good()) << "error in csv read in";
  getline(relativeData, value, ','); // initial comma
  std::vector<Time> measTimes;
  std::vector<ValueType> measValues;
  Eigen::VectorXd m(7);
  double val;
  SE3 pose;
  int nMeas = 401;
  int counterMeas = 0;
  while (relativeData.good() && counterMeas< nMeas) {
    getline (relativeData, value, ',');
    timestamp = atof(value.c_str());
    measTimes.push_back(timestamp);

    for (int i=0; i<7; ++i){
      getline (relativeData, value, ',');
      val = atof(value.c_str());
      m[i] = val;
    }
    pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
    measValues.push_back(pose);
    counterMeas++;
  }
  relativeData.close();

  // fit curve to the coefficients (only time spacing is important here, but constructor demands coefficients as well)
  // create coefficients, one more than the number of measurements
  vector<Time> coefTimes;
  vector<ValueType> coefValues;
  for(int i=0; i <= measTimes.size(); i++) {
    coefTimes.push_back(i);
    coefValues.push_back(SE3(SO3(1, 0, 0, 0),(SE3::Position(0,0,0))));
  }
  SlerpSE3Curve curve;
  curve.fitCurve(coefTimes, coefValues);

  // noise models
  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas((gtsam::Vector(DIM) << 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)); //noise model for GPS
  gtsam::noiseModel::Diagonal::shared_ptr priorNoise = gtsam::noiseModel::Diagonal::
      Sigmas((gtsam::Vector(DIM) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)); //noise model for prior

  // factor graph
  gtsam::NonlinearFactorGraph graph;

  // prior factor (an SE3AbsolutePoseFactor at the good first coefficient)
  SE3 prior(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));

  Expression<ChartValue<ValueType> > predictedPrior(convertToChartValue<ValueType>,
                                                    curve.getEvalExpression2(coefTimes[0]));

  ExpressionFactor<ChartValue<ValueType> > f(priorNoise,
                                             ChartValue<ValueType>(prior),
                                             predictedPrior);
  graph.add(f);

//  for (int i=0; i < measValues.size(); ++i) {
//
//
//    Expression<ChartValue<ValueType> > relative(::relativeMeasurementExpression, curve.getEvalExpression2(measTimes[i]), curve.getEvalExpression(measTimes[i]+1));
//    ExpressionFactor<ChartValue<ValueType> > factor(measNoiseModel,ChartValue<ValueType>(measValues[i]), relative);
//
//    graph.add(factor);
//  }
//
//  // Populate GTSAM values
//  KeyCoefficientTime *rval0, *rval1;
//  Values initials, expected;
//  for(int i=0; i< coefValues.size(); i++) {
//    curve.getCoefficientsAt(coefTimes[i], &rval0, &rval1);
//    Key key;
//    // todo use another method for getting the keys
//    // here a if is necessary since t=maxtime the coef is in rval1
//    if (i == coefValues.size() - 1) {
//      key = rval1->key;
//    } else {
//      key = rval0->key;
//    }
//    initials.insert(key,ValueType(SO3(1,SO3::Vector3(0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0)));
//    //expected.insert(key,valuesCoef[i]);
//  }
//
//  // optimize the trajectory
//  gtsam::LevenbergMarquardtParams params;
//  params.setVerbosity("ERROR");
//  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials, params).optimize();
//
//  for(int i=0; i<coefValues.size(); i++) {
//    curve.getCoefficientsAt(coefTimes[i], &rval0, &rval1);
//    Key key;
//    // todo use another method for getting the keys
//    // here a if is necessary since t=maxtime the coef is in rval1
//    if (i == coefValues.size() - 1) {
//      key = rval1->key;
//    } else {
//      key = rval0->key;
//    }
//    Eigen::Vector3d val = result.at<ValueType>(key).getPosition();
//    cout << val[0] << " " << val[1] << " " << val[2] << endl;
//  }

//  ASSERT_TRUE(expected.equals(result, 0.5));
}
