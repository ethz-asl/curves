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
#include <curves/SE3CoefficientImplementation.hpp>
#include "test_Helpers.hpp"

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
#include <math.h>

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
  bool recordCsv = false;

  // parameters for LM optimizer
  gtsam::LevenbergMarquardtParams params;
  params.setRelativeErrorTol(1e-5);
  params.setlambdaUpperBound(1e5);
  params.setlambdaLowerBound(1e-5);
  params.setAbsoluteErrorTol(-1e5);
  params.setErrorTol(1);
  params.setMaxIterations(50);
  params.setVerbosity("ERROR");

  // filename/location for storing resultfiles
  std::string filename = "/home/johnny/projects/ctsm/src/curves/curves/test/dump/optimized_coefficients_MITb";

  typedef SlerpSE3Curve::ValueType ValueType;
  typedef SlerpSE3Curve::EvaluatorTypePtr EvaluatorTypePtr;

  /// **************** \\\
  ///   DATA READ-IN   \\\
  /// **************** \\\

  // read the absolute measurements (initials) from file
  std::vector<curves::Time> measTimesVertex, measTimesEdge, measTimesLoopA, measTimesLoopB, expectedTimesVertex;
  std::vector<Eigen::VectorXd> outValues;
  std::vector<ValueType> measValuesVertex, measValuesEdge, measValuesLoop, expectedValuesVertex;
  CurvesTestHelpers::loadTimeVectorCSV("MITb_vertex_initials.csv", &measTimesVertex, &outValues);
  SE3 pose;
  for (int i=0; i<outValues.size();++i){
    pose = SE3(SO3(outValues[i][4], outValues[i][5], outValues[i][6], outValues[i][7]),
               (SE3::Position() << outValues[i][1], outValues[i][2], outValues[i][3]).finished());
    measValuesVertex.push_back(pose);
  }
  // read the relative measurements (odometry) from file
  CurvesTestHelpers::loadTimeVectorCSV("MITb_edges_measurements.csv", &measTimesEdge, &outValues);
  for (int i=0; i<outValues.size();++i){
    pose = SE3(SO3(outValues[i][4], outValues[i][5], outValues[i][6], outValues[i][7]),
               (SE3::Position() << outValues[i][1], outValues[i][2], outValues[i][3]).finished());
    measValuesEdge.push_back(pose);
  }
  // read the loop measurements from file
  CurvesTestHelpers::loadTimeTimeVectorCSV("MITb_loop_closure.csv", &measTimesLoopA, &measTimesLoopB, &outValues);
  for (int i=0; i<outValues.size();++i){
    pose = SE3(SO3(outValues[i][4], outValues[i][5], outValues[i][6], outValues[i][7]),
               (SE3::Position() << outValues[i][1], outValues[i][2], outValues[i][3]).finished());
    measValuesLoop.push_back(pose);
  }
  // read the expected values from file
  CurvesTestHelpers::loadTimeVectorCSV("MITb_expected_result.csv", &expectedTimesVertex, &outValues);
  for (int i=0; i<outValues.size();++i){
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
  for(int i=0; i < measTimesVertex.size(); i++) {
    coefTimes.push_back(measTimesVertex[i]);
    coefValues.push_back(measValuesVertex[i]);
    expectedValues.push_back(expectedValuesVertex[i]);
  }
  SlerpSE3Curve curve;
  curve.fitCurve(coefTimes, coefValues);

  // Populate GTSAM values
  KeyCoefficientTime *rval0, *rval1;
  Values initials, expected;
  for(int i=0; i< coefValues.size(); i++) {
    curve.getCoefficientsAt(coefTimes[i], &rval0, &rval1);
    Key key;
    // todo use another method for getting the keys
    // here a if is necessary since t=maxtime the coef is in rval1
    if (i == coefValues.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    initials.insert(key,coefValues[i]);
    expected.insert(key,expectedValues[i]);
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
  Expression<ChartValue<ValueType> > predictedPrior(convertToChartValue<ValueType>,
                                                    curve.getEvalExpression2(coefTimes[0]));
  ExpressionFactor<ChartValue<ValueType> > f(priorNoise,
                                             ChartValue<ValueType>(prior),
                                             predictedPrior);
  graph.add(f);

  // odometry
  for (int i=0; i < measValuesEdge.size(); ++i) {
    Expression<ValueType> TA(curve.getEvalExpression2(measTimesEdge[i]));
    Expression<ValueType> TB(curve.getEvalExpression2(measTimesEdge[i]+1));
    Expression<ValueType> iTA(curves::inverseTransformation, TA);
    Expression<ValueType> relative(curves::composeTransformations, iTA, TB);
    Expression<ChartValue<ValueType> >predicted(gtsam::convertToChartValue<ValueType>, relative);
    ExpressionFactor<ChartValue<ValueType> > factor(measNoiseModel,ChartValue<ValueType>(measValuesEdge[i]), predicted);
    graph.add(factor);
  }

  // loop closures
  for (int i = 0; i < measValuesLoop.size(); ++i) {
    Expression<ValueType> TA(curve.getEvalExpression2(measTimesLoopA[i]));
    Expression<ValueType> TB(curve.getEvalExpression2(measTimesLoopB[i]));
    Expression<ValueType> iTA(curves::inverseTransformation, TA);
    Expression<ValueType> relative(curves::composeTransformations, iTA, TB);
    Expression<ChartValue<ValueType> >predicted(gtsam::convertToChartValue<ValueType>, relative);
    ExpressionFactor<ChartValue<ValueType> > factor(loopNoiseModel,ChartValue<ValueType>(measValuesLoop[i]), predicted);
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
      for(int i2=0; i2<coefValues.size(); ++i2) {
        curve.getCoefficientsAt(coefTimes[i2], &rval0, &rval1);
        Key key;
        // todo use another method for getting the keys
        // here a if is necessary since t=maxtime the coef is in rval1
        if (i2 == coefValues.size() - 1) {
          key = rval1->key;
        } else {
          key = rval0->key;
        }
        Eigen::Vector3d val = initials.at<ValueType>(key).getPosition();
        SO3 val2 = initials.at<ValueType>(key).getRotation();
        resultFile << measTimesVertex[i] << ", "
                   << val2.w() << ", " << val2.x() << ", " << val2.y() << ", " << val2.z()<< ", "
                   << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
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
    CHECK(initials.equals(expected, 1)) << "results don't match the expected";
}
