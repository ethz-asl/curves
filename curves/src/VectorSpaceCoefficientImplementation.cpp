#include <curves/VectorSpaceCoefficientImplementation.hpp>
#include <iostream>
#include <glog/logging.h>

namespace curves {

VectorSpaceCoefficientImplementation::VectorSpaceCoefficientImplementation(unsigned dimension) :
    dimension_(dimension){};

VectorSpaceCoefficientImplementation::~VectorSpaceCoefficientImplementation() {}
  
bool VectorSpaceCoefficientImplementation::equals(const Eigen::VectorXd& thisCoeff, 
                                                  const Eigen::VectorXd& otherCoeff, 
                                                  double tol) const {
  return (thisCoeff - otherCoeff).array().abs().maxCoeff() < tol;
}

Eigen::VectorXd* VectorSpaceCoefficientImplementation::makeUniqueInPlace(Eigen::VectorXd* thisCoeff ) const { 
  CHECK_NOTNULL(thisCoeff);
  return thisCoeff;
}

void VectorSpaceCoefficientImplementation::makeUnique(const Eigen::VectorXd& thisCoeff,
                                                      Eigen::VectorXd* outUniqueCoeff) const {
  CHECK_NOTNULL(outUniqueCoeff);
  *outUniqueCoeff = thisCoeff;
}
 
void VectorSpaceCoefficientImplementation::print(const Eigen::VectorXd& thisCoeff, 
                                                 const std::string& str) const {
  std::cout << str << " " << thisCoeff.transpose();
}

size_t VectorSpaceCoefficientImplementation::dim() const {
  return dimension_;
}

size_t VectorSpaceCoefficientImplementation::ambientDim() const {
  return dimension_;
}

void VectorSpaceCoefficientImplementation::retract(const Eigen::VectorXd& thisCoeff, 
                                                   const Eigen::VectorXd& delta, 
                                                   Eigen::VectorXd* outIncrementedCoeff) const {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outIncrementedCoeff);
  *outIncrementedCoeff = thisCoeff + delta;
}

void VectorSpaceCoefficientImplementation::localCoordinates(const Eigen::VectorXd& thisCoeff, 
                                                            const Eigen::VectorXd& otherCoeff, 
                                                            Eigen::VectorXd* outLocalCoordinates) const {
  // \todo PTF Add check for dimension.
  CHECK_NOTNULL(outLocalCoordinates);
  *outLocalCoordinates = otherCoeff - thisCoeff;
}


} // namespace curves
