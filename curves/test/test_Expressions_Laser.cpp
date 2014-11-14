/*
 * @file test_Expressions_Laser.cpp
 * @date Nov 15, 2014
 * @author Renaud Dube, Abel Gawel
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

#define DIM 6

using namespace curves;
using namespace gtsam;
using namespace std;

typedef typename SlerpSE3Curve::ValueType ValueType;
typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;

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


// Wrapper to enable numerical differentiation on SlerpSE3Evaluator::evaluate
class ExpressionValueWrapper {
 private:
  Expression<ValueType> expression_;
  const Values& initialCoefficients_;

 public:
  ExpressionValueWrapper(const Expression<ValueType>& expression,
                         const Values& initialCoefficients) : expression_(expression),
                         initialCoefficients_(initialCoefficients){

  }
  ~ExpressionValueWrapper() {}

  ChartValue<ValueType> evaluate(const ChartValue<ValueType>& x, int arg) {
    std::set<Key>::const_iterator it = expression_.keys().begin();
    if (arg) {++it;}
    Key key = *it;

    Values values(initialCoefficients_);
    values.update(key,x.value());

    std::vector<size_t> dimensions;
    dimensions.push_back(DIM);
    dimensions.push_back(DIM);
    static const int Dim = traits::dimension<ValueType>::value;
    VerticalBlockMatrix Ab(dimensions, Dim);
    FastVector<Key> keys = boost::assign::list_of(*(expression_.keys().begin()))(*(++expression_.keys().begin()));
    JacobianMap actualMap(keys,Ab);

    ValueType result = expression_.value(values, actualMap);

    return convertToChartValue<ValueType>(result);
  }
};


TEST(CurvesTestSuite, spiralLaser) {
  SlerpSE3Curve curve;
  const double t[] = {0, 90, 180, 181, 225, 270, 360};
  //  const double evalTime = 5;

  // setup two SE3 objects
  // ValueType poseA(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(0,0,0));
  ValueType poseA(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(1,0,0));
  ValueType poseB(SO3(0.7071067811865476,SO3::Vector3(0.7071067811865475,0.0,0)),SE3::Position(2,0,0));
  ValueType poseC(SO3(0.0,SO3::Vector3(1,0.0,0)),SE3::Position(3,0,0));
  ValueType poseD(SO3(0.0,SO3::Vector3(-1,0.0,0)),SE3::Position(3,0,0));
  ValueType poseE(SO3(0.38268343236508967,SO3::Vector3(-0.923879532511287,0,0)),SE3::Position(3.5,0,0));
  ValueType poseF(SO3(0.7071067811865475,SO3::Vector3(-0.7071067811865475,0.0,0)),SE3::Position(4,0,0));
  ValueType poseG(SO3(1,SO3::Vector3(0,0,0)),SE3::Position(5,0,0));

  std::vector<Time> times(t,t+7);
  std::vector<ValueType> values;
  values.push_back(poseA);
  values.push_back(poseB);
  values.push_back(poseC);
  values.push_back(poseD);
  values.push_back(poseE);
  values.push_back(poseF);
  values.push_back(poseG);
  curve.fitCurve(times, values);



  // fake laserdata:
  double const distance = 5;
  double const angleIncrement = 1; // deg
  double const numMeas = 360;

  std::vector<Time> timesMeas;
  for (int i=0; i<numMeas; ++i) {
    timesMeas.push_back(i);
  }

  // Populate GTSAM values
  KeyCoefficientTime *rval0, *rval1;
  gtsam::Values gtsamValues;
  for(int i=0; i< values.size(); i++) {
    curve.getCoefficientsAt(times[i], &rval0, &rval1);
    Key key;
    if (i == values.size() - 1) {
      key = rval1->key;
    } else {
      key = rval0->key;
    }
    gtsamValues.insert(key,values[i]);
  }

  const Eigen::Vector3d pos(0,distance,0);
  std::vector<Eigen::Vector3d> pointCloud;
  std::ofstream resultFile;
  resultFile.open("laser_spiral.csv");
  std::cout << "pc: " <<std::endl;
  for (int i=0; i<numMeas; ++i) {
    Expression<ValueType> expression1 = curve.getEvalExpression2(timesMeas[i]);
    SE3 val = expression1.value(gtsamValues);
//    std::cout << val.getPosition()[0] << " " <<val.getPosition()[1]  << " "<< val.getPosition()[2] <<std::endl;
//    std::cout << val.getRotation()[0] << " " <<val.getRotation()[1]  << " "<< val.getRotation()[2] << " "<< val.getRotation()[3] <<std::endl;
    if (i >175 && i< 185) {

    std::cout << val.getRotation() <<std::endl;
    }
    pointCloud.push_back(val.transform(pos));
    resultFile<< i << ", " << pointCloud[i][0] << ", " <<pointCloud[i][1]  << ", "<< pointCloud[i][2]<< std::endl;

//    std::cout << pointCloud[i][0] << " " <<pointCloud[i][1]  << " "<< pointCloud[i][2] <<std::endl;
  }
  resultFile.close();

}


