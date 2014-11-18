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

// test a greater example of using more factors, creating gtsam factor graph and optimizing it
TEST(CurvesTestSuite, test_MITb) {

  typedef SlerpSE3Curve::ValueType ValueType;
  typedef SlerpSE3Curve::EvaluatorTypePtr EvaluatorTypePtr;

  // read the absolute measurements from file
  curves::Time timestampVertex;
  std::string valueVertex;
  std::ifstream absoluteData ("MITb_vertex_initials.csv");
  CHECK(absoluteData.good()) << "error in csv read in";
  getline(absoluteData, valueVertex, ','); // initial comma
  std::vector<Time> measTimesVertex;
  std::vector<ValueType> measValuesVertex;
  Eigen::VectorXd m(13);
  double val;
  SE3 pose;
  int nMeasVertex = 808;
  int counterMeas = 0;
  while (absoluteData.good() && counterMeas< nMeasVertex) {
    getline (absoluteData, valueVertex, ',');
    timestampVertex = atof(valueVertex.c_str());
    measTimesVertex.push_back(timestampVertex);

    for (int i=0; i<7; ++i){
      getline (absoluteData, valueVertex, ',');
      val = atof(valueVertex.c_str());
      m[i] = val;
    }
    pose = SE3(SO3(m[3], m[4], m[5], m[6]),(SE3::Position() << m[0],m[1],m[2]).finished());
    measValuesVertex.push_back(pose);
    counterMeas++;
  }
  absoluteData.close();

  // read the relative measurements from file
  curves::Time timestamp;
  std::string value;
  std::ifstream relativeData ("MITb_edges_measurements.csv");
  CHECK(relativeData.good()) << "error in csv read in";
  getline(relativeData, value, ','); // initial comma
  Eigen::Matrix<double,6,1> measCov;
  std::vector<Eigen::Matrix<double,6,1> > measCovariances;
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
    measCov << m[7], m[8], m[9], m[10], m[11], m[12];
    measCovariances.push_back(measCov);
    counterMeas++;
  }
  relativeData.close();


  // read the loop closure measurements from file
  std::ifstream loopData ("MITb_loop_closure.csv");
  CHECK(loopData.good()) << "error in csv read in";
  getline(loopData, value, ','); // initial comma
  std::vector<Time> loopTimesA, loopTimesB;
  curves::Time timestampA, timestampB;
  std::vector<ValueType> loopValues;
  Eigen::Matrix<double,6,1> loopCov;
  std::vector<Eigen::Matrix<double,6,1> > loopCovariances;
  int nLoop = 20;
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
    loopCov << m[7], m[8], m[9], m[10], m[11], m[12];
    loopCovariances.push_back(loopCov);

    counterMeas++;
  }
  loopData.close();

  // fit curve to the coefficients (only time spacing is important here, but constructor demands coefficients as well)
  // create coefficients, one more than the number of measurements
  vector<Time> coefTimes;
  vector<ValueType> coefValues;
  for(int i=0; i < measTimesVertex.size(); i++) {
    coefTimes.push_back(measTimesVertex[i]);
    coefValues.push_back(measValuesVertex[i]);
    std::cout << measTimesVertex[i] <<std::endl;
    //    std::cout << "Position: " << i << " "<< measValuesVertex[i].getPosition() << std::endl;
    //    std::cout << "Rotation: " << i << " " << measValuesVertex[i].getRotation() <<std::endl;
    //    coefValues.push_back(SE3(SO3(1, 0, 0, 0),(SE3::Position(0,0,0))));
  }


  SlerpSE3Curve curve;
  curve.fitCurve(coefTimes, coefValues);

  // noise models
//  gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
//      Sigmas((gtsam::Vector(DIM) << 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)); //noise model for GPS
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

  // \todo give reasonable prior
  graph.add(f);

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
    //        initials.insert(key,ValueType(SO3(1,SO3::Vector3(0.0,0.0,0.0)),SE3::Position(0.0,0.0,0.0)));
    //expected.insert(key,valuesCoef[i]);
    //        initials.print("initials:");

    //    std::cout << "diff: " << relativeMeasurementExpressionX(coefValues[i],coefValues[i+1]).getPosition() <<std::endl;
  }
  for (int i=0; i < measValues.size(); ++i) {


    //    Expression<ChartValue<ValueType> > relative(::relativeMeasurementExpressionX, curve.getEvalExpression2(measTimes[i]), curve.getEvalExpression2(measTimes[i]+1));
    gtsam::noiseModel::Diagonal::shared_ptr measNoiseModel = gtsam::noiseModel::Diagonal::
              Sigmas((gtsam::Vector(DIM) << (measCovariances[i])[5], (measCovariances[i])[3], 0, 0, 0, (measCovariances[i])[1])); //noise model for GPS

    Expression<ValueType> TA(curve.getEvalExpression2(measTimes[i]));
    Expression<ValueType> TB(curve.getEvalExpression2(measTimes[i]+1));
    Expression<ValueType> iTA(curves::inverseTransformation, TA);
    Expression<ValueType> relative(curves::composeTransformations, iTA, TB);

    Expression<ChartValue<ValueType> >predicted(gtsam::convertToChartValue<ValueType>, relative);
    ExpressionFactor<ChartValue<ValueType> > factor(measNoiseModel,ChartValue<ValueType>(measValues[i]), predicted);
    //    std::cout << "i: " << i << " "<< factor.unwhitenedError(initials).transpose() <<std::endl;
    graph.add(factor);
  }

  for (int i = 0; i<loopValues.size(); ++i) {
//  for (int i = 10; i<16; ++i) {

//    std::cout << "covariances: " <<(loopCovariances[i])[0] << " " << (loopCovariances[i])[3] << " " <<(loopCovariances[i])[5] << std::endl;
    gtsam::noiseModel::Diagonal::shared_ptr loopNoiseModel = gtsam::noiseModel::Diagonal::
          Sigmas((gtsam::Vector(DIM) << 5*(loopCovariances[i])[5], (loopCovariances[i])[3], 0, 0, 0, (loopCovariances[i])[1])); //noise model for GPS
//    gtsam::Matrix6 covariance;
//    covariance << (loopCovariances[i])[0],(loopCovariances[i])[1],0,0,0,0,
//                (-loopCovariances[i])[1], (loopCovariances[i])[3], 0, 0, 0, 0,
//                0, 0, 0, 0, 0, 0,
//                0, 0, 0, 0, 0, 0,
//                0, 0, 0, 0, 0, 0,
//                0, 0, 0, 0, 0, (loopCovariances[i])[5];
//    gtsam::noiseModel::Gaussian::shared_ptr loopNoiseModel = gtsam::noiseModel::Gaussian::Covariance
//        (covariance); //noise model for GPS
    Expression<ValueType> TA(curve.getEvalExpression2(loopTimesA[i]));
    Expression<ValueType> TB(curve.getEvalExpression2(loopTimesB[i]));
    Expression<ValueType> iTA(curves::inverseTransformation, TA);
    Expression<ValueType> relative(curves::composeTransformations, iTA, TB);

    Expression<ChartValue<ValueType> >predicted(gtsam::convertToChartValue<ValueType>, relative);
    ExpressionFactor<ChartValue<ValueType> > factor(loopNoiseModel,ChartValue<ValueType>(loopValues[i]), predicted);
    graph.add(factor);
    factor.print("loop factor");
  }
  //initials.print("initials:");
  //graph.print("graph");

  //  std::cout << "c6: " << coefValues[6].getPosition() <<std::endl;
  //    std::cout << "c7: " << coefValues[7].getPosition() <<std::endl;
  //    std::cout << "diff 6 - 7: " << relativeMeasurementExpressionX(coefValues[6], coefValues[7]).getPosition() <<std::endl;



  std::cout << __LINE__ << " " << __FILE__ << std::endl;

  //
  // optimize the trajectory
  gtsam::LevenbergMarquardtParams params;
  params.setVerbosity("ERROR");
  //  params.setVerbosity("LINEAR");
  gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials, params).optimize();

  std::ofstream resultFile;
  resultFile.open("optimized_coefficients_MITb.csv");
  for(int i=0; i<coefValues.size(); i++) {
    curve.getCoefficientsAt(coefTimes[i], &rval0, &rval1);
    Key key;
    // todo use another method for getting the keys
    // here a if is necessary since t=maxtime the coef is in rval1
    if (i == coefValues.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    Eigen::Vector3d val = result.at<ValueType>(key).getPosition();
    resultFile << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
  }
  //    cout << val[0] << " " << val[1] << " " << val[2] << endl;
  //  ASSERT_TRUE(expected.equals(result, 0.5));
}
