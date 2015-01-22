/*
 * @file test_SE3_Sinus_Circle.cpp
 * @date Nov 05, 2014
 * @author Abel Gawel, Renaud Dube
 *
 *  Note : this test should eventually be transfered to trajectory maintainer
 */

#include <gtest/gtest.h>
#include <curves/SlerpSE3Curve.hpp>
#include "gtsam/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include "gtsam/nonlinear/ExpressionFactor.h"
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
typedef SlerpSE3Curve::CoefficientIter CoefficientIter;

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

// test a greater example of using more factors, creating gtsam factor graph and optimizing it
TEST(CurvesTestSuite, testSE3AbsolutePoseFactor_3333_SOxR3_GPS) {

  // load the true coefficients' values for comparison
  double val;
  curves::Time timestamp;
  std::string value;
  std::ifstream coeffData ( "3D_sinus_circle_coefficients.csv" );
  CHECK(coeffData.good()) << "error in csv read in";
  getline (coeffData, value, ','); // initial comma
  std::vector<Time> times;
  std::vector<ValueType> valuesCoef, initialsCoef;
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
    initialsCoef.push_back(SE3(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0)));
    counterCoef++;
  }
  coeffData.close();
  //cout << "ncoef " << valuesCoef.size();

  // fit curve to the coefficients (only time spacing is important here, but constructor demands coefficients as well)
  SlerpSE3Curve curve;
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(times, initialsCoef, &outKeys);

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

  // noise models
  Vector6 gpsNoise;
  gpsNoise << 0.05, 0.05, 0.05, 0.05, 0.05, 0.05;
  Vector6 priNoise;
  priNoise << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

  noiseModel::Diagonal::shared_ptr GPSNoise = noiseModel::Diagonal::Sigmas(gpsNoise);
  noiseModel::Diagonal::shared_ptr priorNoise = noiseModel::Diagonal::Sigmas(priNoise);

  // factor graph
  gtsam::NonlinearFactorGraph graph;

  // prior factor (an SE3AbsolutePoseFactor at the good first coefficient)
  SE3 prior(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));

  Expression<ValueType> predictedPrior = curve.getValueExpression(times[0]);

  ExpressionFactor<ValueType> f(priorNoise, prior, predictedPrior);
  graph.add(f);

  for (size_t i=0; i < valuesGPS.size(); ++i) {
    Expression<ValueType> predicted = curve.getValueExpression(timesGPS[i]);

    ExpressionFactor<ValueType> f(GPSNoise, valuesGPS[i], predicted);
    graph.add(f);
  }

  // Populate GTSAM values
  CoefficientIter rval0, rval1;
  Values initials, expected;
  curve.initializeGTSAMValues(&initials);
  for(size_t i=0; i< outKeys.size(); i++) {
    expected.insert(outKeys[i],valuesCoef[i]);
  }

  // optimize the trajectory
  gtsam::LevenbergMarquardtParams params;
  // params.setVerbosity("ERROR");
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials, params).optimize();

  ASSERT_TRUE(expected.equals(result, 0.5));
}
