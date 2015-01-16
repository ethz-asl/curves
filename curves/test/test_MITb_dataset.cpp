/*
 * @file test_MITb_dataset.cpp
 * @date Nov 05, 2014
 * @author Abel Gawel, Renaud Dube
 *
 *  Note : this test should eventually be transfered to trajectory maintainer
 *  Dataset source: www.lucacarlone.com/index.php/resources/datasets
 *
 *  This test performs pose graph SLAM on the MITb dataset, using GTSAM Expressions
 *  and curve library
 *
 *  See DEFINITIONS section for adjusting parameters
 */

#include <gtest/gtest.h>
#include <curves/SlerpSE3Curve.hpp>
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

#define DIM 6
#define OUTPUT_RESULTS_AS_CSV

using namespace curves;
using namespace std;
using namespace gtsam;

typedef SlerpSE3Curve::ValueType ValueType;

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;

TEST(CurvesTestSuite, test_MITb) {

  /// **************** \\\
  ///    DEFINITIONS   \\\
  /// **************** \\\

  gtsam::Vector6 loopNoise;
  loopNoise << 0.1, 0.1, 0.001, 0.2*M_PI/180, 0.2*M_PI/180, 1*M_PI/180;
  gtsam::Vector6 measNoise;
  measNoise << 0.1, 0.1, 0.001, 0.1*M_PI/180, 0.1*M_PI/180, 0.5*M_PI/180;
  gtsam::Vector6 priNoise;
  priNoise << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  bool recordCsv = true;

  // parameters for LM optimizer
  gtsam::LevenbergMarquardtParams params;
  params.setRelativeErrorTol(1e-5);
  params.setlambdaUpperBound(1e5);
  params.setlambdaLowerBound(1e-5);
  params.setAbsoluteErrorTol(-1e5);
  params.setErrorTol(1);
  params.setMaxIterations(50);
  params.setVerbosity("SILENT");

  // filename/location for storing resultfiles
  std::string filename = "/home/renaud/projects/ctsm/src/curves/curves/test/dump/optimized_coefficients_MITb";

  typedef SlerpSE3Curve::ValueType ValueType;

  /// **************** \\\
  ///   DATA READ-IN   \\\
  /// **************** \\\

  // read the absolute measurements (initials) from file
  std::vector<curves::Time> measTimesVertex, measTimesEdge, measTimesLoopA, measTimesLoopB, expectedTimesVertex;
  std::vector<Eigen::VectorXd> outValues;
  std::vector<ValueType> measValuesVertex, measValuesEdge, measValuesLoop, expectedValuesVertex;
  CurvesTestHelpers::loadTimeVectorCSV("MITb_vertex_initials.csv", &measTimesVertex, &outValues);
  SE3 pose;
  for (size_t i=0; i<outValues.size();++i){
    pose = SE3(SO3(outValues[i][4], outValues[i][5], outValues[i][6], outValues[i][7]),
               (SE3::Position() << outValues[i][1], outValues[i][2], outValues[i][3]).finished());
    measValuesVertex.push_back(pose);
  }
  // read the relative measurements (odometry) from file
  CurvesTestHelpers::loadTimeVectorCSV("MITb_edges_measurements.csv", &measTimesEdge, &outValues);
  for (size_t i=0; i<outValues.size();++i){
    pose = SE3(SO3(outValues[i][4], outValues[i][5], outValues[i][6], outValues[i][7]),
               (SE3::Position() << outValues[i][1], outValues[i][2], outValues[i][3]).finished());
    measValuesEdge.push_back(pose);
  }
  // read the loop measurements from file
  CurvesTestHelpers::loadTimeTimeVectorCSV("MITb_loop_closure.csv", &measTimesLoopA, &measTimesLoopB, &outValues);
  for (size_t i=0; i<outValues.size();++i){
    pose = SE3(SO3(outValues[i][4], outValues[i][5], outValues[i][6], outValues[i][7]),
               (SE3::Position() << outValues[i][1], outValues[i][2], outValues[i][3]).finished());
    measValuesLoop.push_back(pose);
  }
  // read the expected values from file
  CurvesTestHelpers::loadTimeVectorCSV("MITb_expected_result.csv", &expectedTimesVertex, &outValues);
  for (size_t i=0; i<outValues.size();++i){
    pose = SE3(SO3(outValues[i][0], outValues[i][1], outValues[i][2], outValues[i][3]),
               (SE3::Position() << outValues[i][4], outValues[i][5], outValues[i][6]).finished());
    expectedValuesVertex.push_back(pose);
  }

  /// **************** \\\
  ///   CREATE CURVE   \\\
  /// **************** \\\

  // fit curve to the coefficients (only time spacing is important here, but constructor demands coefficients as well)
  // create coefficients, one more than the number of measurements
  vector<Time> coefTimes;
  vector<ValueType> coefValues, expectedValues;
  for(size_t i=0; i < measTimesVertex.size(); i++) {
    coefTimes.push_back(measTimesVertex[i]);
    coefValues.push_back(measValuesVertex[i]);
    expectedValues.push_back(expectedValuesVertex[i]);
  }
  SlerpSE3Curve curve;
  std::vector<gtsam::Key> outKeys;
  curve.fitCurve(coefTimes, coefValues, &outKeys);

  // Populate GTSAM values
  Values initials, expected;
  curve.initializeGTSAMValues(&initials);
  for(size_t i=0; i< outKeys.size(); i++) {
    expected.insert(outKeys[i],expectedValues[i]);
  }

  /// **************** \\\
  ///   NOISE MODELS   \\\
  /// **************** \\\

  //noise model for prior
  gtsam::noiseModel::Diagonal::shared_ptr priorNoise = gtsam::noiseModel::Diagonal::
      Sigmas(priNoise);
  // noise model for odometry measurements
  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(measNoise);
  // noise model for loop closure ("odometry"-like) measurements
  gtsam::noiseModel::Diagonal::shared_ptr loopNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(loopNoise);

  /// **************** \\\
  ///      FACTORS     \\\
  /// **************** \\\

  // factor graph
  gtsam::NonlinearFactorGraph graph;

  // prior
  SE3 prior(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  Expression<ValueType> predictedPrior = curve.getValueExpression2(coefTimes[0]);
  ExpressionFactor<ValueType> f(priorNoise, prior, predictedPrior);
  graph.add(f);

  // odometry
  for (size_t i=0; i < measValuesEdge.size(); ++i) {
    Expression<ValueType> TA(curve.getValueExpression2(measTimesEdge[i]));
    Expression<ValueType> TB(curve.getValueExpression2(measTimesEdge[i]+1));
    Expression<ValueType> predicted = kindr::minimal::invertAndCompose(TA,TB);
    ExpressionFactor<ValueType> factor(measNoiseModel,measValuesEdge[i], predicted);
    graph.add(factor);
  }

  // loop closures
  for (size_t i = 0; i < measValuesLoop.size(); ++i) {
    Expression<ValueType> TA(curve.getValueExpression2(measTimesLoopA[i]));
    Expression<ValueType> TB(curve.getValueExpression2(measTimesLoopB[i]));
    Expression<ValueType> predicted = kindr::minimal::invertAndCompose(TA,TB);
    ExpressionFactor<ValueType> factor(loopNoiseModel,measValuesLoop[i], predicted);
    graph.add(factor);
  }

  // iterate to incrementally narrow down noise of loop closures
  for (int i=0; i<12; ++i) {

    /// **************** \\\
    ///  RECORD RESULTS  \\\
    /// **************** \\\

    if(recordCsv) {
      std::stringstream ss;
      ss << filename << i << ".csv";
      std::ofstream resultFile;
      resultFile.open(ss.str().c_str());

      for (size_t i=0; i< outKeys.size(); i++) {
        Eigen::Vector3d pos = initials.at<ValueType>(outKeys[i]).getPosition();
        SO3 rot = initials.at<ValueType>(outKeys[i]).getRotation();
        resultFile << measTimesVertex[i] << ", "
            << rot.w() << ", " << rot.x() << ", " << rot.y() << ", " << rot.z()<< ", "
            << pos[0] << ", " << pos[1] << ", " << pos[2] << std::endl;
      }
    }

    /// **************** \\\
    ///   OPTIMIZATION   \\\
    /// **************** \\\

    // catch last iteration (only for writing results, if activated)
    if(i<11) {
      // narrow down noise model for loop closures
      gtsam::noiseModel::Diagonal::shared_ptr loopNoiseModelX = gtsam::noiseModel::Diagonal::
          Sigmas(loopNoise*(101-(i*10)));
      // update noise model of loop closures
      *loopNoiseModel = *loopNoiseModelX;
      // optimize graph
      initials = gtsam::LevenbergMarquardtOptimizer(graph, initials, params).optimize();
    }
  }
  ASSERT_TRUE(initials.equals(expected, 5)) << "results don't match the expected";
}

ExpressionFactor<ValueType>
getRelativeMeasurementFactor(const SlerpSE3Curve& curve,
                             Time timeA, Time timeB,
                             ValueType measurement,
                             noiseModel::Diagonal::shared_ptr noiseModel) {
  Expression<ValueType> TA(curve.getValueExpression2(timeA));
  Expression<ValueType> TB(curve.getValueExpression2(timeB));
  Expression<ValueType> predicted = kindr::minimal::invertAndCompose(TA,TB);
  return ExpressionFactor<ValueType>(noiseModel, measurement, predicted);
}

// This test uses the ISAM2 algorithm to optimize over the MITb dataset.
// Coefficients are added incrementally to the ISAM Bayes tree.
// Optionally, results at each iteration can be saved to a .csv file for visualization.
TEST(CurvesTestSuite, test_MITb_ISAM2) {

  typedef SlerpSE3Curve::ValueType ValueType;

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
  SE3 pose;
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
    pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
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
    pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
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
    pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
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
    pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
    expectedValues.push_back(pose);
  }
  expectedData.close();

  // Fit curve to the coefficients
  SlerpSE3Curve curve;
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
  gtsam::Vector6 measNoise;
  measNoise << 0.1, 0.1, 0.001, 0.1*M_PI/180, 0.1*M_PI/180, 0.5*M_PI/180;
  gtsam::Vector6 loopNoise;
  loopNoise << 0.2, 0.2, 0.002, 0.2*M_PI/180, 0.2*M_PI/180, 1*M_PI/180;
  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(measNoise);
  gtsam::noiseModel::Diagonal::shared_ptr loopNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(loopNoise);

  // Create a Factor Graph and gtsam structures to hold initial and expected values
  NonlinearFactorGraph graph;
  Values gtsamInitial, gtsamExpected;

  // Loop over all coefficients, adding the observations to iSAM incrementally
  for (size_t i = 0; i < coefficientTimes.size(); ++i) {

    // Get the desired key which will be pushed to ISAM2 as an update
    // Key keyToAdd = mapTimes.at(coefficientTimes[i]);

    // Build odometry measurement factors and push them to the graph
    for (int z = 0; z < measTimes.size(); z++) {
      // Before, we used :
      // curve.getCoefficientsAt(measTimes[z], &rval0, &rval1);
      // if ( rval0->key == keyToAdd) { add factor }
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
    for (int z = 0; z < loopTimesA.size(); z++) {
      //      curve.getCoefficientsAt(loopTimesA[z], &rval0, &rval1);
      //      if ( rval0->key == keyToAdd) {
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
      SE3 prior(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
      Expression<ValueType> predictedPrior = curve.getValueExpression2(coefficientTimes[0]);
      ExpressionFactor<ValueType> f(priorNoise, prior, predictedPrior);
      graph.push_back(f);
    } else {
      // Create gtsam initial values to be pushed to ISAM2 and expected results.
      FastVector<Key> keys = graph.keys();
      FastVector<Key>::const_iterator it = keys.begin();
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
      ss << "/home/renaud/MATLAB/MITb/optimized_coefficients_MITb" << i << ".csv";

      std::ofstream resultFile;
      resultFile.open(ss.str().c_str());
      for(int i2=0; i2<=i; i2++) {
        Eigen::Vector3d pos = currentEstimate.at<ValueType>(outKeys[i2]).getPosition();
        Eigen::Vector4d rot = currentEstimate.at<ValueType>(outKeys[i2]).getRotation().vector();
        resultFile << coefficientTimes[i2] << ", " << pos[0] << ", " << pos[1] << ", " << pos[2] << ", " <<
            rot[0] << ", " << rot[1] << ", " << rot[2] << ", " << rot[3] << std::endl;
      }
#endif

      // Check the solution when all coefficients have been added to the problem
      if (i == coefficientTimes.size()-1) {
        ASSERT_TRUE(currentEstimate.equals(gtsamExpected, 0.5));
      }

      // Clear the factor graph and values for the next iteration
      graph.resize(0);
      gtsamInitial.clear();
    }
  }
}
