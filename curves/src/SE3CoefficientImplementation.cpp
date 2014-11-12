#include <curves/SE3CoefficientImplementation.hpp>
#include <iostream>
#include <glog/logging.h>
#include "kindr/minimal/quat-transformation.h"

namespace curves {

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

void SE3CoefficientImplementation::makeValue(const SE3& pose, Eigen::VectorXd *outValue) const {
  CHECK_NOTNULL(outValue);
  CHECK_EQ(7,outValue->size());
  (*outValue) << pose.getPosition(),pose.getRotation().vector();
}

void SE3CoefficientImplementation::retract(const Eigen::VectorXd& thisCoeff,
                                           const Eigen::VectorXd& delta,
                                           Eigen::VectorXd* outIncrementedCoeff) const {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outIncrementedCoeff);

//  // SO x R3 retract
//  SO3 thisSO3(SO3::Vector4(thisCoeff.segment<4>(3)));
//  SE3::Position thisSE3Position = thisCoeff.head<3>();
//  SO3 updatedRotation = SO3(delta.tail<3>().eval())*thisSO3;
//  SE3::Position updatedTranslation = SE3::Position(delta.head<3>().eval()) + thisSE3Position;
//  (*outIncrementedCoeff) << updatedTranslation, updatedRotation.vector();

  // SE3 retract
  // the position is stored in the first 3 dimenstions, and the quaternion is in the next 4 or the coeff vector
  SE3 thisSE3(SO3(thisCoeff[3], thisCoeff.segment<3>(4)),thisCoeff.head<3>());
//  SE3 thisSE3(SO3(SO3::Vector4(thisCoeff.segment<4>(3))),thisCoeff.head<3>());
  // the SE3 constructor with a 6D vector is the exponential map
  SE3 updated = SE3(delta.head<6>().eval())*thisSE3;
  (*outIncrementedCoeff) << updated.getPosition(), updated.getRotation().vector();
}

// chart-friendly overload of retract
void SE3CoefficientImplementation::retract(const SE3& thisSE3,
                                           const Eigen::Matrix<double,6,1>& delta,
                                           SE3* outIncrementedSE3) {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outIncrementedSE3);

  // the position is stored in the first 3 dimenstions, and the quaternion is in the next 4 or the coeff vector
    // the SE3 constructor with a 6D vector is the exponential map
  SE3 updated = SE3(delta)*thisSE3;
  (*outIncrementedSE3) = updated;
}

void SE3CoefficientImplementation::localCoordinates(const Eigen::VectorXd& thisCoeff,
                                                    const Eigen::VectorXd& otherCoeff,
                                                    Eigen::VectorXd* outLocalCoordinates) const {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outLocalCoordinates);

//  // SO x R3 local coordinates
//  SO3 thisSO3(SO3::Vector4(thisCoeff.segment<4>(3)));
//  SO3 otherSO3(SO3::Vector4(otherCoeff.segment<4>(3)));
//  SE3::Position thisPosition(thisCoeff.head<3>());
//  SE3::Position otherPosition(otherCoeff.head<3>());
//  SO3 deltaRotation = otherSO3 * thisSO3.inverted();
//  AngleAxis deltaRotationAA(deltaRotation);
//  SE3::Position deltaPosition = otherPosition - thisPosition;
//  (*outLocalCoordinates) << deltaPosition, deltaRotationAA.axis()*deltaRotationAA.angle();

  // SE3 local coordinates
  SE3 thisSE3(SO3(thisCoeff[3], thisCoeff.segment<3>(4)),thisCoeff.head<3>());
  SE3 otherSE3(SO3(otherCoeff[3], otherCoeff.segment<3>(4)),otherCoeff.head<3>());
  // local coordinates are defined to be on the left side of the transformation,
  // ie other = [delta]^ A
  SE3 delta = otherSE3*thisSE3.inverted();
  (*outLocalCoordinates) = delta.log();
}

// chart-friendly overload of localCoordinates
void SE3CoefficientImplementation::localCoordinates(const SE3& thisSE3,
                                                    const SE3& otherSE3,
                                                    Eigen::Matrix<double,6,1>* outLocalCoordinates) {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outLocalCoordinates);
  // SE3 local coordinates
  // local coordinates are defined to be on the left side of the transformation,
  // ie other = [delta]^ A
  SE3 delta = otherSE3*thisSE3.inverted();
  (*outLocalCoordinates) = delta.log();
}

} // namespace curves




