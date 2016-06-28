/*
 * @file test_MITb_dataset.cpp
 * @date Nov 25, 2015
 * @author Renaud Dubé, Abel Gawel
 *
 *  Dataset source: www.lucacarlone.com/index.php/resources/datasets
 *
 *  This test performs pose graph SLAM on the MITb dataset, using GTSAM Expressions
 *  and curve library
 *
 *  See DEFINITIONS section for adjusting parameters
 */

#include <gtest/gtest.h>
#include <curves/SlerpSE2Curve.hpp>
#include <curves/Pose2_Expressions.hpp>
#include "test_Helpers.hpp"

#include "gtsam/nonlinear/Expression.h"
#include <gtsam/nonlinear/NonlinearFactorGraph.h>
#include <gtsam/nonlinear/LevenbergMarquardtOptimizer.h>
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include "gtsam/nonlinear/ExpressionFactor.h"
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/base/Value.h>
#include <gtsam/base/GenericValue.h>
#include <gtsam/nonlinear/ISAM2.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#define DIM 3
#define OUTPUT_RESULTS_AS_CSV

using namespace curves;
using namespace std;
using namespace gtsam;

typedef SlerpSE2Curve::ValueType ValueType;

typedef Pose2 SE2;
typedef Rot2 SO2;

double convertQuatToZAngle(double q0, double q1, double q2, double q3) {
  return atan2(2*(q0*q3+q1*q2),1-2*(q2*q2+q3*q3));
}

//
ExpressionFactor<ValueType>
getRelativeMeasurementFactor(const SlerpSE2Curve& curve,
                             Time timeA, Time timeB,
                             ValueType measurement,
                             noiseModel::Diagonal::shared_ptr noiseModel) {
  Expression<ValueType> TA(curve.getValueExpression(timeA));
  Expression<ValueType> TB(curve.getValueExpression(timeB));
  Expression<ValueType> predicted = compose(inverse(TA),TB);
  return ExpressionFactor<ValueType>(noiseModel, measurement, predicted);
}

// This test uses the ISAM2 algorithm to optimize over the MITb dataset.
// Coefficients are added incrementally to the ISAM Bayes tree.
// Optionally, results at each iteration can be saved to a .csv file for visualization.
TEST(CurvesTestSuite, test_MITb_ISAM2_SE2) {

  typedef SlerpSE2Curve::ValueType ValueType;

  // Read the initial estimates from file
  curves::Time timestampVertex;
  std::string valueVertex;
  std::ifstream absoluteData ("MITb_vertex_initials.csv");
  CHECK(absoluteData.good()) << "error in csv read in";
  getline(absoluteData, valueVertex, ','); // initial comma
  std::vector<Time> coefficientTimes;
  std::vector<ValueType> initialValues;
  Eigen::VectorXd m(13);
  double val;
  SE2 pose;
  int nMeasVertex = 808;
  int counterMeas = 0;
  while (absoluteData.good() && counterMeas< nMeasVertex) {
    getline (absoluteData, valueVertex, ',');
    timestampVertex = atof(valueVertex.c_str());
    coefficientTimes.push_back(timestampVertex);
    for (int i=0; i<7; ++i){
      getline (absoluteData, valueVertex, ',');
      val = atof(valueVertex.c_str());
      m[i] = val;
    }
    pose = SE2(m[0], m[1], convertQuatToZAngle(m[3], m[4], m[5], m[6]));
    initialValues.push_back(pose);
    counterMeas++;
  }
  absoluteData.close();

  // Read the relative measurements from file
  curves::Time timestamp;
  std::string value;
  std::ifstream relativeData ("MITb_edges_measurements.csv");
  CHECK(relativeData.good()) << "error in csv read in";
  getline(relativeData, value, ','); // initial comma
  std::vector<Time> measTimes;
  std::vector<ValueType> measValues;
  int nMeas = 807;
  counterMeas = 0;
  while (relativeData.good() && counterMeas< nMeas) {
    getline (relativeData, value, ',');
    timestamp = atof(value.c_str());
    measTimes.push_back(timestamp);

    for (int i=0; i<13; ++i){
      getline (relativeData, value, ',');
      val = atof(value.c_str());
      m[i] = val;
    }
    pose = SE2(m[0], m[1], convertQuatToZAngle(m[3], m[4], m[5], m[6]));
    measValues.push_back(pose);
    counterMeas++;
  }
  relativeData.close();

  // Read the loop closure measurements from file
  std::ifstream loopData ("MITb_loop_closure.csv");
  CHECK(loopData.good()) << "error in csv read in";
  getline(loopData, value, ','); // initial comma
  std::vector<Time> loopTimesA, loopTimesB;
  curves::Time timestampA, timestampB;
  std::vector<ValueType> loopValues;
  int nLoop = 15;
  counterMeas = 0;
  while (loopData.good() && counterMeas< nLoop) {
    getline (loopData, value, ',');
    timestampA = atof(value.c_str());
    loopTimesA.push_back(timestampA);
    getline (loopData, value, ',');
    timestampB = atof(value.c_str());
    loopTimesB.push_back(timestampB);

    for (int i=0; i<13; ++i){
      getline (loopData, value, ',');
      val = atof(value.c_str());
      m[i] = val;
    }
    pose = SE2(m[0], m[1], convertQuatToZAngle(m[3], m[4], m[5], m[6]));
    loopValues.push_back(pose);
    counterMeas++;
  }
  loopData.close();

  // Read the expected result from file
  std::ifstream expectedData ("MITb_expected_ISAM2_result.csv");
  CHECK(expectedData.good()) << "error in csv read in";
  getline(expectedData, value, ','); // initial comma
  std::vector<Time> expectedTimes;
  std::vector<ValueType> expectedValues;
  while (expectedData.good()) {
    getline (loopData, value, ',');
    timestamp = atof(value.c_str());
    expectedTimes.push_back(timestamp);
    for (int i=0; i<7; ++i){
      getline (expectedData, value, ',');
      val = atof(value.c_str());
      m[i] = val;
    }
    pose = SE2(m[0], m[1], convertQuatToZAngle(m[3], m[4], m[5], m[6]));
    expectedValues.push_back(pose);
  }
  expectedData.close();

  // Fit curve to the coefficients
  SlerpSE2Curve curve;
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(coefficientTimes, initialValues, &outKeys);

  // Make maps of Key to ValueType and of Time to Key
  boost::unordered_map<Key, ValueType> mapInitialValues;
  boost::unordered_map<Key, ValueType> mapExpectedValues;
  std::map<Time, Key> mapTimes;
  for(size_t i=0; i< outKeys.size(); i++) {
    mapInitialValues[outKeys[i]] = initialValues[i];
    mapExpectedValues[outKeys[i]] = expectedValues[i];
    mapTimes[coefficientTimes[i]] = outKeys[i];
  }

  // Create an ISAM2 object
  ISAM2Params parameters;
  parameters.relinearizeThreshold = 0.01;
  parameters.relinearizeSkip = 1;
  ISAM2 isam(parameters);

  // Define the noise models
  gtsam::Vector3 measNoise;
  measNoise << 0.1, 0.1, 0.5*M_PI/180;
  gtsam::Vector3 loopNoise;
  loopNoise << 0.2, 0.2, 1*M_PI/180;


  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(measNoise);
  gtsam::noiseModel::Diagonal::shared_ptr loopNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(loopNoise);

  // Create a Factor Graph and gtsam structures to hold initial and expected values
  NonlinearFactorGraph graph;
  Values gtsamInitial, gtsamExpected;

  // Loop over all coefficients, adding the observations to iSAM incrementally
  for (size_t i = 0; i < coefficientTimes.size(); ++i) {

    // Build odometry measurement factors and push them to the graph
    for (size_t z = 0; z < measTimes.size(); z++) {
      // Before, we used :
      bool addFactor = false;
      if (i < coefficientTimes.size() - 1) {
        if ( measTimes[z] >= coefficientTimes[i] && measTimes[z] < coefficientTimes[i+1] ) {
          addFactor = true;
        }
      } else {
        if (measTimes[z] == coefficientTimes[i]) {
          addFactor = true;
        }
      }
      if (addFactor)
        graph.push_back(getRelativeMeasurementFactor(curve, measTimes[z], measTimes[z]+1,
                                                     measValues[z], measNoiseModel));
    }

    // Build loop closure measurement factors and push them to the graph
    for (size_t z = 0; z < loopTimesA.size(); z++) {
      bool addFactor = false;
      if (i < coefficientTimes.size() - 1) {
        if ( loopTimesA[z] >= coefficientTimes[i] && loopTimesA[z] < coefficientTimes[i+1] ) {
          addFactor = true;
        }
      } else {
        if (loopTimesA[z] == coefficientTimes[i]) {
          addFactor = true;
        }
      }
      if (addFactor) {
        graph.push_back(getRelativeMeasurementFactor(curve, loopTimesA[z], loopTimesB[z],
                                                     loopValues[z], loopNoiseModel));
      }
    }

    // If this is the first iteration, add a prior on the first coefficient
    if( i == 0) {
      gtsam::noiseModel::Constrained::shared_ptr priorNoise = gtsam::noiseModel::Constrained::All(DIM);
      SE2 prior(0,0,0);
      Expression<ValueType> predictedPrior = curve.getValueExpression(coefficientTimes[0]);
      ExpressionFactor<ValueType> f(priorNoise, prior, predictedPrior);
      graph.push_back(f);
    } else {
      // Create gtsam initial values to be pushed to ISAM2 and expected results.
      gtsam::KeySet keys = graph.keys();
      gtsam::KeySet::const_iterator it = keys.begin();
      for (; it != keys.end(); ++it) {
        if ( !gtsamExpected.exists(*it) ) {
          gtsamInitial.insert(*it, mapInitialValues.at(*it));
          gtsamExpected.insert(*it, mapExpectedValues.at(*it));
        }
      }

      // Update iSAM with the new factors
      isam.update(graph, gtsamInitial);
      // Each call to iSAM2 update(*) performs one iteration of the iterative nonlinear solver.
      // If accuracy is desired at the expense of time, update(*) can be called additional times
      // to perform multiple optimizer iterations every step.
      isam.update();

      Values currentEstimate = isam.calculateEstimate();

#ifdef OUTPUT_RESULTS_AS_CSV
      std::stringstream ss;
      ss << "/home/renaud/MATLAB/MITb/SE2/optimized_coefficients_MITb" << i << ".csv";

      std::ofstream resultFile;
      resultFile.open(ss.str().c_str());
      for(size_t i2=0; i2<=i; i2++) {
        SE2 val = currentEstimate.at<ValueType>(outKeys[i2]);
        resultFile << coefficientTimes[i2] << ", " << val.x() << ", " << val.y() << ", " << val.theta() << std::endl;
      }
#endif

      // Check the solution when all coefficients have been added to the problem
      if (i == coefficientTimes.size()-1) {
        ASSERT_TRUE(currentEstimate.equals(gtsamExpected, 5));
      }

      // Clear the factor graph and values for the next iteration
      graph.resize(0);
      gtsamInitial.clear();
    }
  }
}
