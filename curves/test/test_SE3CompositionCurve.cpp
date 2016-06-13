/*
 * @file test_SE3CombinationCurve.hpp
 * @date Feb 06, 2015
 * @author Abel Gawel, Renaud Dubé, Mike Bosse
 */

#include <gtest/gtest.h>
#include <curves/SlerpSE3Curve.hpp>
#include <curves/CubicHermiteSE3Curve.hpp>
#include <curves/SE3CompositionCurve.hpp>
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

#include <kindr/minimal/testing-gtsam.h>

#include "test_Helpers.hpp"

#include <eigen-checks/gtest.h>

using namespace curves;
using namespace gtsam;
using namespace std;

template <typename CurveType_, boost::int64_t MinSamplingPeriod_>
class CurveWrapper {
 public:
  typedef CurveType_ CurveType;

 private:
  CurveType curve_;

 public:
  CurveWrapper() {
    curve_.setMinSamplingPeriod(MinSamplingPeriod_);
  }

  CurveType getCurve() {
    return curve_;
  }
};

template <typename CurveWrapperType_>
class SinusCircleTestSuites : public ::testing::Test {
 public:

  typedef typename CurveWrapperType_::CurveType CurveType;

  CurveType createCurve() {
    CurveWrapperType_ curveWrapper;
    return curveWrapper.getCurve();
  }
};

typedef ::testing::Types<
    CurveWrapper<SE3CompositionCurve<SlerpSE3Curve, SlerpSE3Curve>, 1000000000>,
    CurveWrapper<SE3CompositionCurve<SlerpSE3Curve, SlerpSE3Curve>, 500000000>
> TestTypes;

TYPED_TEST_CASE(SinusCircleTestSuites, TestTypes);

typedef typename SlerpSE3Curve::ValueType ValueType;
typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;

template <class Curve>
ExpressionFactor<ValueType>
getRelativeMeasurementFactor(const Curve& curve,
                             const Time timeA, const Time timeB,
                             const ValueType& measurement,
                             gtsam::noiseModel::Diagonal::shared_ptr noiseModel,
                             Values values) {
  Expression<ValueType> TA(curve.getValueExpression(timeA));
  Expression<ValueType> TB(curve.getValueExpression(timeB));
  Expression<ValueType> iTA(kindr::minimal::inverse(TA));
  Expression<ValueType> relative(kindr::minimal::compose(iTA, TB));

  gtsam::testExpressionJacobians(relative, values, 1e-9, 1e-4);

  return ExpressionFactor<ValueType>(noiseModel,measurement,relative);
}

// todo move to another file / place
double multiplyVectors(Eigen::Vector3d a, Eigen::Vector3d b,
                       OptionalJacobian<1,3> H) {

  if(H)
    *H = a.transpose();

  return a.transpose() * b;
}

TYPED_TEST(SinusCircleTestSuites, testSE3CompositionCurve_SinusCircle) {

  // Load the true coefficients
  double val;
  curves::Time timestamp;
  std::string value;
  std::ifstream coeffData ( "3D_sinus_circle_coefficients_composition.csv" );
  CHECK(coeffData.good()) << "error in csv read in";
  getline (coeffData, value, ','); // initial comma
  std::vector<Time> times;
  std::vector<ValueType> valuesCoef, initialsCoef;
  Eigen::VectorXd m(7);
  SE3 pose;
  int counterCoef = 0;
  while (coeffData.good()) {
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
  }
  coeffData.close();

  // Load the matches
  std::ifstream matchesData ( "3D_sinus_circle_matches.csv" );
  CHECK(matchesData.good()) << "error in csv read in";
  getline (matchesData, value, ','); // initial comma
  Eigen::VectorXd m2(8);
  std::vector<Eigen::VectorXd> matches;
  while (matchesData.good()) {
    for (int i=0; i<8; ++i){
      getline (matchesData, value, ',');
      val = atof(value.c_str());
      m2[i] = val;
    }
    matches.push_back(m2);
  }
  matchesData.close();

  // Make the true curve
  SlerpSE3Curve trueCurve;
  trueCurve.extend(times, valuesCoef);

  // Build noisy odometry data
  const int nBaseCoefToConsider = valuesCoef.size();
  std::vector<ValueType> noisyOdom;
  for (int i = 0; i < nBaseCoefToConsider - 1; ++i) {
    ValueType TA(valuesCoef[i]);
    ValueType TB(valuesCoef[i+1]);
    ValueType T_A_B = TA.inverted() * TB;

    ValueType noise;
    noise.setRandom(0.0003, 0.0003);
    ValueType bias(SE3::Position(0.0003,0,0),SE3::Rotation(1,0,0,0));
    noisyOdom.push_back(bias*T_A_B*noise);
  }

  std::vector<ValueType> noisyInitials;
  noisyInitials.push_back(valuesCoef[0]);
  for (int i = 0; i < nBaseCoefToConsider - 1; ++i) {
    SE3 integrated(noisyInitials[i] * noisyOdom[i]);
    integrated = SE3(integrated.getRotation().normalize(), integrated.getPosition());
    noisyInitials.push_back(integrated);
  }

  // Make the composition curve
  typename TypeParam::CurveType curve = this->createCurve();
  for (size_t i=0; i < nBaseCoefToConsider; ++i) {
    std::vector<Time> timeToAdd;
    std::vector<ValueType> valueToAdd;
    timeToAdd.push_back(times[i]);
    valueToAdd.push_back(noisyInitials[i]);
    curve.extend(timeToAdd, valueToAdd);
  }

  // noise models
  Vector6 measNoise;
  measNoise << 0.005, 0.005, 0.005, 0.008727, 0.008727, 0.008727;
  Vector6 priNoise;
  priNoise << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
  Vector6 loopClosureNoise;
  loopClosureNoise << 0.001, 0.001, 0.001, 0.001, 0.001, 0.001;
  noiseModel::Diagonal::shared_ptr measNoiseModel = noiseModel::Diagonal::Sigmas(measNoise);
  noiseModel::Diagonal::shared_ptr priNoiseModel = noiseModel::Diagonal::Sigmas(priNoise);
  noiseModel::Diagonal::shared_ptr loopClosureNoiseModel = noiseModel::Diagonal::Sigmas(loopClosureNoise);
  Eigen::Matrix<double,1,1> laserNoise;
  laserNoise << 0.1;
  gtsam::noiseModel::Diagonal::shared_ptr laserNoiseModel = gtsam::noiseModel::Diagonal::
      Sigmas(laserNoise);

  // Optimization loop
  for (int i = 0; i < 1; ++i) {

    // Populate GTSAM values
    Values initials;
    curve.initializeGTSAMValues(&initials);

    // factor graph
    gtsam::NonlinearFactorGraph graph;

    // prior factor (an SE3AbsolutePoseFactor at the good first coefficient)
    Expression<ValueType> predictedPrior = curve.getValueExpression(times[0]);
    ExpressionFactor<ValueType> f(priNoiseModel, valuesCoef[0], predictedPrior);
    graph.add(f);

    // Add the odometry factors
    for (size_t i = 0; i < nBaseCoefToConsider - 1; ++i) {
      graph.add(getRelativeMeasurementFactor(curve, times[i], times[i+1],
                                             noisyOdom[i],
                                             measNoiseModel, initials));
    }

    // Add the match factors
    for (size_t i = 0; i < matches.size(); ++i) {

      Time meanTimeA = matches[i][0];
      Time meanTimeB = matches[i][1];

      Eigen::Vector3d mu_W_SA(matches[i][2],matches[i][3],matches[i][4]);
      mu_W_SA += Eigen::Vector3d::Random()/100;
      Eigen::Vector3d mu_W_SB(matches[i][2],matches[i][3],matches[i][4]);
      mu_W_SA += Eigen::Vector3d::Random()/100;

      const ValueType true_T_W_A_const = trueCurve.evaluate(meanTimeA);   // here take from real values?
      const ValueType true_T_W_B_const = trueCurve.evaluate(meanTimeB);

      const Eigen::Vector3d mu_A_SA  = true_T_W_A_const.inverted().transform(mu_W_SA);
      const Eigen::Vector3d mu_B_SB  = true_T_W_B_const.inverted().transform(mu_W_SB);

      Expression<Eigen::Vector3d> E_mu_A_SA(mu_A_SA);
      Expression<Eigen::Vector3d> E_mu_B_SB(mu_B_SB);

      Expression<SE3> T_W_B(curve.getValueExpression(meanTimeB));
      Expression<Eigen::Vector3d> pointTransformedB = kindr::minimal::transform(T_W_B, E_mu_B_SB);

      Expression<SE3> T_W_A(curve.getValueExpression(meanTimeA));
      Expression<Eigen::Vector3d> pointTransformedA = kindr::minimal::transform(T_W_A, E_mu_A_SA);

      Expression<Eigen::Vector3d> substracted = kindr::minimal::vectorDifference(pointTransformedA, pointTransformedB);

      Eigen::Vector3d norm(matches[i][5],matches[i][6],matches[i][7]);

      Expression<double> error(boost::bind(&multiplyVectors,norm,_1,_2),substracted);
      ExpressionFactor<double> matchFactor(laserNoiseModel,(double (0)), error);
      graph.push_back(matchFactor);
    }

    ///////////////
    // Save the graph
    std::filebuf fb;
    fb.open ("/tmp/test.dot",std::ios::out);
    std::ostream os(&fb);
    graph.saveGraph(os);

    ///////////////
    // Save the graph connectivity
    FastSet<Key> keysInGraph = graph.keys();

    std::vector<curves::Time> factorIds;
    std::vector<Eigen::VectorXd> factorsKeysToSave;

    factorIds.push_back(0);

    Eigen::VectorXd keysVector(keysInGraph.size());
    int counter = 0;
    for(FastSet<Key>::const_iterator it = keysInGraph.begin(); it != keysInGraph.end(); ++it) {
      keysVector[counter] = *it;
      counter++;
    }
    factorsKeysToSave.push_back(keysVector);

    for(gtsam::NonlinearFactorGraph::const_iterator it = graph.begin(); it != graph.end(); ++it) {
      gtsam::FastVector<Key> keysInFactor = (*it)->keys();
      keysVector.setZero();
      for (gtsam::FastVector<Key>::const_iterator it2 = keysInFactor.begin(); it2 < keysInFactor.end(); ++it2) {
        keysVector[std::distance(keysInGraph.begin(), keysInGraph.find(*it2))] = 1;

      }
      factorIds.push_back(factorIds[factorIds.size() -1] + 1);
      factorsKeysToSave.push_back(keysVector);

    }
    CurvesTestHelpers::writeTimeVectorCSV("/tmp/factors.csv", factorIds, factorsKeysToSave);

    // Perform optimization
    gtsam::LevenbergMarquardtParams params;
    gtsam::Values result = gtsam::LevenbergMarquardtOptimizer(graph, initials, params).optimize();

    curve.updateFromGTSAMValues(result);

    ///////////////
    // Save residuals before & after optimization
    std::vector<Time> ids;
    std::vector<Eigen::VectorXd> errorsBefore, errorsAfter;
    Eigen::VectorXd error(1);
    Time counter1 = 0;
    for(gtsam::NonlinearFactorGraph::const_iterator it = graph.begin(); it != graph.end(); ++it) {
      ids.push_back(counter);
      counter1++;
      error << (*it)->error(initials);
      errorsBefore.push_back(error);
      error << (*it)->error(result);
      errorsAfter.push_back(error);
    }
    CurvesTestHelpers::writeTimeVectorCSV("/tmp/circle_sinus/residualsBeforeOptimizationComposition.csv", ids, errorsBefore);
    CurvesTestHelpers::writeTimeVectorCSV("/tmp/circle_sinus/residualsAfterOptimizationComposition.csv", ids, errorsAfter);

    ///////////////
    // Save the curve values

    std::vector<Eigen::VectorXd> expected, initial, optimized;
    std::vector<Time> timesToSave;
    for (unsigned int z = 0; z < nBaseCoefToConsider; ++z){
      timesToSave.push_back(times[z]);

      SE3 val;

      Eigen::VectorXd v1(7), v2(7), v3(7);

      val = valuesCoef[z];
      v1 << val.getPosition().x(), val.getPosition().y(), val.getPosition().z(),
          val.getRotation().w(), val.getRotation().x(), val.getRotation().y(), val.getRotation().z();
      expected.push_back(v1);

      val = noisyInitials[z];
      v2 << val.getPosition().x(), val.getPosition().y(), val.getPosition().z(),
          val.getRotation().w(), val.getRotation().x(), val.getRotation().y(), val.getRotation().z();
      initial.push_back(v2);

      val = curve.evaluate(times[z]);
      v3 << val.getPosition().x(), val.getPosition().y(), val.getPosition().z(),
          val.getRotation().w(), val.getRotation().x(), val.getRotation().y(), val.getRotation().z();
      optimized.push_back(v3);
    }

    CurvesTestHelpers::writeTimeVectorCSV("/tmp/circle_sinus/expected.csv", timesToSave, expected);
    CurvesTestHelpers::writeTimeVectorCSV("/tmp/circle_sinus/initial.csv", timesToSave, initial);
    CurvesTestHelpers::writeTimeVectorCSV("/tmp/circle_sinus/optimized.csv", timesToSave, optimized);
  }

}
