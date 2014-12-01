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
#include "gtsam/nonlinear/Expression.h"
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


// test a greater example of using more factors, creating gtsam factor graph and optimizing it
TEST(CurvesTestSuite, testSE3AbsolutePoseFactor_3333_SOxR3_GPS) {

  typedef SlerpSE3Curve::ValueType ValueType;
  typedef SlerpSE3Curve::EvaluatorTypePtr EvaluatorTypePtr;

  // load the true coefficients' values for comparison
  double val;
  curves::Time timestamp;
  std::string value;
  std::ifstream coeffData ( "3D_sinus_circle_coefficients.csv" );
  CHECK(coeffData.good()) << "error in csv read in";
  getline (coeffData, value, ','); // initial comma
  std::vector<Time> times;
  std::vector<ValueType> valuesCoef;
  Eigen::VectorXd m(7);
  SE3 pose;
  int nCoef = 401;
  int counterCoef = 0;
  while (coeffData.good() && counterCoef < nCoef) {
    getline (coeffData, value, ',');
    timestamp = atof(value.c_str());
    times.push_back(timestamp);
    for (int i=0; i<7; ++i){
      getline (coeffData, value, ',');
      val = atof(value.c_str());
      m[i] = val;
    }
    pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
    valuesCoef.push_back(pose);
    counterCoef++;
  }
  coeffData.close();
  //cout << "ncoef " << valuesCoef.size();

  // fit curve to the coefficients (only time spacing is important here, but constructor demands coefficients as well)
  SlerpSE3Curve curve;
  curve.fitCurve(times, valuesCoef);

  // read the GPS measurements from file
  std::ifstream GPSData ("3D_sinus_circle_GPS_clean.csv");
  CHECK(GPSData.good()) << "error in csv read in";
  getline (GPSData, value, ','); // initial comma
  std::vector<Time> timesGPS;
  std::vector<ValueType> valuesGPS;
  bool reachedLastCoef = false;
  while (GPSData.good() && !reachedLastCoef) {
    getline (GPSData, value, ',');
    timestamp = atof(value.c_str());

    if (timestamp <= times[nCoef-1]) {
      timesGPS.push_back(timestamp);
      for (int i=0; i<7; ++i){
        getline (GPSData, value, ',');
        val = atof(value.c_str());
        m[i] = val;
      }
      pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
      valuesGPS.push_back(pose);
    } else {
      reachedLastCoef = true;
    }

  }
  GPSData.close();
  //cout << "valuesGPS " << valuesGPS.size();

  // noise models
  gtsam::noiseModel::Diagonal::shared_ptr GPSNoise = gtsam::noiseModel::Diagonal::
      Sigmas((gtsam::Vector(DIM) << 0.05, 0.05, 0.05, 0.05, 0.05, 0.05)); //noise model for GPS
  gtsam::noiseModel::Diagonal::shared_ptr priorNoise = gtsam::noiseModel::Diagonal::
      Sigmas((gtsam::Vector(DIM) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)); //noise model for prior

  // factor graph
  gtsam::NonlinearFactorGraph graph;

  // prior factor (an SE3AbsolutePoseFactor at the good first coefficient)
  SE3 prior(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));

  Expression<ChartValue<ValueType> > predictedPrior(convertToChartValue<ValueType>,
                                                    curve.getEvalExpression2(times[0]));

  ExpressionFactor<ChartValue<ValueType> > f(priorNoise,
                                             ChartValue<ValueType>(prior),
                                             predictedPrior);
  graph.add(f);

  for (int i=0; i < valuesGPS.size(); ++i) {
    Expression<ChartValue<ValueType> > predicted(convertToChartValue<ValueType>,
                                                 curve.getEvalExpression2(timesGPS[i]));

    ExpressionFactor<ChartValue<ValueType> > f(GPSNoise,
                                               ChartValue<ValueType>(valuesGPS[i]),
                                               predicted);
    graph.add(f);
  }

  // Populate GTSAM values
  KeyCoefficientTime *rval0, *rval1;
  Values initials, expected;
  for(int i=0; i< valuesCoef.size(); i++) {
    curve.getCoefficientsAt(times[i], &rval0, &rval1);
    Key key;
    // todo use another method for getting the keys
    // here a if is necessary since t=maxtime the coef is in rval1
    if (i == valuesCoef.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    initials.insert(key,ValueType(SO3(1,SO3::Vector3(0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0)));
    expected.insert(key,valuesCoef[i]);
  }

  // optimize the trajectory
  gtsam::LevenbergMarquardtParams params;
  // params.setVerbosity("ERROR");
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials, params).optimize();

  ASSERT_TRUE(expected.equals(result, 0.5));
}
