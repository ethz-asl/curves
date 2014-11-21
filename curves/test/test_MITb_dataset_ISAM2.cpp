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
#include <gtsam/nonlinear/GaussNewtonOptimizer.h>
#include "gtsam_unstable/nonlinear/ExpressionFactor.h"
#include <gtsam/base/numericalDerivative.h>
#include <gtsam/nonlinear/ISAM2.h>

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

#define DIM 6
// #define OUTPUT_RESULTS_AS_CSV

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


ExpressionFactor<ChartValue<ValueType> >
getRelativeMeasurementFactor(const SlerpSE3Curve& curve,
                             Time timeA, Time timeB,
                             ValueType measurement,
                             gtsam::noiseModel::Diagonal::shared_ptr noiseModel) {
  Expression<ValueType> TA(curve.getEvalExpression2(timeA));
  Expression<ValueType> TB(curve.getEvalExpression2(timeB));
  Expression<ValueType> iTA(curves::inverseTransformation, TA);
  Expression<ValueType> relative(curves::composeTransformations, iTA, TB);

  Expression<ChartValue<ValueType> >predicted(gtsam::convertToChartValue<ValueType>, relative);
  return ExpressionFactor<ChartValue<ValueType> >(noiseModel,ChartValue<ValueType>(measurement), predicted);
}

// This test uses the ISAM2 algorithm to optimize over the MITb dataset.
// Coefficients are added incrementally to the ISAM Bayes tree.
// Optionally, results at each iteration can be saved to a .csv file for visualization.
TEST(CurvesTestSuite, test_MITb_ISAM2) {

  typedef SlerpSE3Curve::ValueType ValueType;
  typedef SlerpSE3Curve::EvaluatorTypePtr EvaluatorTypePtr;

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
  curve.fitCurve(coefficientTimes, initialValues);

  // Make maps of Key to ValueType and of Time to Key
  boost::unordered_map<Key, ValueType> mapInitialValues;
  boost::unordered_map<Key, ValueType> mapExpectedValues;
  std::map<Time, Key> mapTimes;
  KeyCoefficientTime *rval0, *rval1;
  for(int i=0; i< initialValues.size(); i++) {
    curve.getCoefficientsAt(coefficientTimes[i], &rval0, &rval1);
    Key key;
    // todo use another method for getting the keys
    // here a if is necessary since t=maxtime the coef is in rval1
    if (i == initialValues.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    mapInitialValues[key] = initialValues[i];
    mapExpectedValues[key] = expectedValues[i];
    mapTimes[coefficientTimes[i]] = key;
  }

  // Create an ISAM2 object
  ISAM2Params parameters;
  parameters.relinearizeThreshold = 0.01;
  parameters.relinearizeSkip = 1;
  ISAM2 isam(parameters);

  // Define the noise models
  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas((gtsam::Vector(DIM) << 0.1, 0.1, 0.001, 0.1*M_PI/180, 0.1*M_PI/180, 0.5*M_PI/180));

  gtsam::noiseModel::Diagonal::shared_ptr loopNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas((gtsam::Vector(DIM) << 0.2, 0.2, 0.002, 0.2*M_PI/180, 0.2*M_PI/180, 1*M_PI/180));

  // Create a Factor Graph and gtsam structures to hold initial and expected values
  NonlinearFactorGraph graph;
  Values gtsamInitial, gtsamExpected;

  // Loop over all coefficients, adding the observations to iSAM incrementally
  for (size_t i = 0; i < coefficientTimes.size(); ++i) {

    // Get the desired key which will be pushed to ISAM2 as an update
    Key keyToAdd = mapTimes.at(coefficientTimes[i]);

    // Build odometry measurement factors and push them to the graph
    for (int z = 0; z < measTimes.size(); z++) {
      curve.getCoefficientsAt(measTimes[z], &rval0, &rval1);
      if ( rval0->key == keyToAdd) {
        graph.push_back(getRelativeMeasurementFactor(curve, measTimes[z], measTimes[z]+1,
                                                     measValues[z], measNoiseModel));
      }
    }

    // Build loop closure measurement factors and push them to the graph
    for (int z = 0; z < loopTimesA.size(); z++) {
      curve.getCoefficientsAt(loopTimesA[z], &rval0, &rval1);
      if ( rval0->key == keyToAdd) {
        graph.push_back(getRelativeMeasurementFactor(curve, loopTimesA[z], loopTimesB[z],
                                                     loopValues[z], loopNoiseModel));
      }
    }

    // If this is the first iteration, add a prior on the first coefficient
    if( i == 0) {
      gtsam::noiseModel::Constrained::shared_ptr priorNoise = gtsam::noiseModel::Constrained::All(DIM);
      SE3 prior(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
      Expression<ChartValue<ValueType> > predictedPrior(convertToChartValue<ValueType>,
                                                        curve.getEvalExpression2(coefficientTimes[0]));
      ExpressionFactor<ChartValue<ValueType> > f(priorNoise,
                                                 ChartValue<ValueType>(prior),
                                                 predictedPrior);
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
        curve.getCoefficientsAt(coefficientTimes[i2], &rval0, &rval1);
        Key key;
        // todo use another method for getting the keys
        // here a if is necessary since t=maxtime the coef is in rval1
        if (i2 == coefficientTimes.size() - 1) {
          key = rval1->key;
        } else {
          key = rval0->key;
        }
        Eigen::Vector3d pos = currentEstimate.at<ValueType>(key).getPosition();
        Eigen::Vector4d rot = currentEstimate.at<ValueType>(key).getRotation().vector();
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
