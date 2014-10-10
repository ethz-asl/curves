#include <curves/SE3CoefficientImplementation.hpp>
#include <iostream>
#include <glog/logging.h>
#include "kindr/minimal/quat-transformation.h"

namespace curves {

SE3CoefficientImplementation::SE3CoefficientImplementation(unsigned dimension) :
    dimension_(dimension){};

SE3CoefficientImplementation::~SE3CoefficientImplementation() {}
  
bool SE3CoefficientImplementation::equals(const Eigen::VectorXd& thisCoeff,
                                                  const Eigen::VectorXd& otherCoeff, 
                                                  double tol) const {
  return (thisCoeff - otherCoeff).array().abs().maxCoeff() < tol;
}

void SE3CoefficientImplementation::makeUniqueInPlace(Eigen::VectorXd* thisCoeff ) const {
  CHECK_NOTNULL(thisCoeff);
}

void SE3CoefficientImplementation::makeUnique(const Eigen::VectorXd& thisCoeff,
                                                      Eigen::VectorXd* outUniqueCoeff) const {
  CHECK_NOTNULL(outUniqueCoeff);
  *outUniqueCoeff = thisCoeff;
}
 
void SE3CoefficientImplementation::print(const Eigen::VectorXd& thisCoeff,
                                                 const std::string& str) const {
  std::cout << str << " " << thisCoeff.transpose();
}

size_t SE3CoefficientImplementation::dim() const {
  return dimension_;
}

size_t SE3CoefficientImplementation::ambientDim() const {
  return dimension_;
}

void SE3CoefficientImplementation::retract(const Eigen::VectorXd& thisCoeff,
                                                   const Eigen::VectorXd& delta, 
                                                   Eigen::VectorXd* outIncrementedCoeff) const {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outIncrementedCoeff);
  *outIncrementedCoeff = thisCoeff + delta;
}

void SE3CoefficientImplementation::localCoordinates(const Eigen::VectorXd& thisCoeff,
                                                            const Eigen::VectorXd& otherCoeff, 
                                                            Eigen::VectorXd* outLocalCoordinates) const {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outLocalCoordinates);
  *outLocalCoordinates = otherCoeff - thisCoeff;
}

} // namespace curves
