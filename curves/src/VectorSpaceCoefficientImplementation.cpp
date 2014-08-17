#include <curves/VectorSpaceCoefficientImplementation.hpp>
#include <iostream>

namespace curves {

VectorSpaceCoefficientImplementation::VectorSpaceCoefficientImplementation(unsigned dimension) :
    dimension_(dimension){};

VectorSpaceCoefficientImplementation::~VectorSpaceCoefficientImplementation() {}
  
bool VectorSpaceCoefficientImplementation::equals(const Eigen::VectorXd& thisCoeff, 
                                                  const Eigen::VectorXd& otherCoeff, 
                                                  double tol) const {
  return (thisCoeff - otherCoeff).array().abs().maxCoeff() < tol;
}

void VectorSpaceCoefficientImplementation::makeUniqueInPlace(Eigen::VectorXd& /* thisCoeff */ ) const{}

void VectorSpaceCoefficientImplementation::makeUnique(const Eigen::VectorXd& thisCoeff,
                                                      Eigen::VectorXd& outUniqueCoeff) const {
  outUniqueCoeff = thisCoeff;
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
                                                   Eigen::VectorXd& outIncrementedCoeff) const {
  // \todo PTF Add check for dimension.
  outIncrementedCoeff = thisCoeff + delta;
}

void VectorSpaceCoefficientImplementation::localCoordinates(const Eigen::VectorXd& thisCoeff, 
                                                            const Eigen::VectorXd& otherCoeff, 
                                                            Eigen::VectorXd& outLocalCoordinates) const {
  // \todo PTF Add check for dimension.
  outLocalCoordinates = otherCoeff - thisCoeff;
}


} // namespace curves
