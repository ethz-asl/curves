#include <curves/SE3CoefficientImplementation.hpp>
#include <iostream>
#include <glog/logging.h>
#include "kindr/minimal/quat-transformation.h"

namespace curves {

typedef kindr::minimal::QuatTransformationTemplate<double> SE3;
typedef SE3::Rotation SO3;
typedef kindr::minimal::AngleAxisTemplate<double> AngleAxis;


SE3CoefficientImplementation::SE3CoefficientImplementation() {}

SE3CoefficientImplementation::~SE3CoefficientImplementation() {}
  
bool SE3CoefficientImplementation::equals(const Eigen::VectorXd& thisCoeff,
                                                  const Eigen::VectorXd& otherCoeff, 
                                                  double tol) const {
  Eigen::VectorXd delta(6);
  localCoordinates(thisCoeff,otherCoeff,&delta);
  return delta.array().abs().maxCoeff() < tol;
}

void SE3CoefficientImplementation::makeUniqueInPlace(Eigen::VectorXd* thisCoeff ) const {
  CHECK_NOTNULL(thisCoeff);
  if ( (*thisCoeff)[3] < 0.0 ) {
    thisCoeff->segment<4>(3) *= -1.0;
  }  
}

void SE3CoefficientImplementation::makeUnique(const Eigen::VectorXd& thisCoeff,
                                                      Eigen::VectorXd* outUniqueCoeff) const {
  CHECK_NOTNULL(outUniqueCoeff);
  *outUniqueCoeff = thisCoeff;
  makeUniqueInPlace(outUniqueCoeff);
}
 
void SE3CoefficientImplementation::print(const Eigen::VectorXd& thisCoeff,
                                                 const std::string& str) const {
  std::cout << str << " " << thisCoeff.transpose();
}

void SE3CoefficientImplementation::makeValue(const Eigen::Matrix4d& matrix, Eigen::VectorXd *outValue) const {
  CHECK_NOTNULL(outValue);
  CHECK_EQ(7,outValue->size());
  SE3 pose(SO3(SO3::RotationMatrix(matrix.topLeftCorner<3,3>())),matrix.topRightCorner<3,1>());
  (*outValue) << pose.getPosition(),pose.getRotation().vector();
}



void SE3CoefficientImplementation::retract(const Eigen::VectorXd& thisCoeff,
                                                   const Eigen::VectorXd& delta, 
                                                   Eigen::VectorXd* outIncrementedCoeff) const {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outIncrementedCoeff);

  SO3 thisSO3(SO3::Vector4(thisCoeff.segment<4>(3)));

    SE3::Position thisSE3Position = thisCoeff.head<3>();

    SO3 updatedRotation = SO3(delta.tail<3>().eval())*thisSO3;
  SE3::Position updatedTranslation = SE3::Position(delta.head<3>().eval()) + thisSE3Position;

  std::cout << __FILE__ <<" : " << __LINE__ <<std::endl;
  // the position is stored in the first 3 dimenstions, and the quaternion is in the next 4 or the coeff vector
//  SE3 thisSE3(SO3(SO3::Vector4(thisCoeff.segment<4>(3))),thisCoeff.head<3>());
  // the SE3 constructor with a 6D vector is the exponential map
//  SE3 updated = SE3(delta.head<6>().eval())*thisSE3;

  (*outIncrementedCoeff) << updatedTranslation, updatedRotation.vector();
}

void SE3CoefficientImplementation::localCoordinates(const Eigen::VectorXd& thisCoeff,
                                                            const Eigen::VectorXd& otherCoeff, 
                                                            Eigen::VectorXd* outLocalCoordinates) const {
  // \todo PTF Add check for dimension.
  std::cout << __FILE__ <<" : " << __LINE__ <<std::endl;
  CHECK_NOTNULL(outLocalCoordinates);


  SO3 thisSO3(SO3::Vector4(thisCoeff.segment<4>(3)));
    SO3 otherSO3(SO3::Vector4(otherCoeff.segment<4>(3)));

    SE3::Position thisPosition(thisCoeff.head<3>());
    SE3::Position otherPosition(otherCoeff.head<3>());

//    SO3 a_R_b = thisSO3.inverted()*otherSO3;

    SO3 deltaRotation = otherSO3 * thisSO3.inverted();
    AngleAxis deltaRotationAA(deltaRotation);
    SE3::Position deltaPosition = otherPosition - thisPosition;
//    delta.setAngle( delta.angle()*alpha_ );
//    SE3 w_T_i(w_R_a*SO3(delta)  , (w_t_a*(1-alpha_)+w_t_b*alpha_).eval());

  std::cout << __FILE__ <<" : " << __LINE__ <<std::endl;
//  SE3 thisSE3(SO3(SO3::Vector4(thisCoeff.segment<4>(3))),thisCoeff.head<3>());
//  SE3 otherSE3(SO3(SO3::Vector4(otherCoeff.segment<4>(3))),otherCoeff.head<3>());

  // local coordinates are defined to be on the left side of the transformation,
  // ie other = [delta]^ A
//  SE3 delta = otherSE3*thisSE3.inverted();

  (*outLocalCoordinates) << deltaPosition, deltaRotationAA.axis()*deltaRotationAA.angle();
}

} // namespace curves
