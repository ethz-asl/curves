#include <curves/CoefficientImplementation.hpp>
#include <iostream>

namespace curves {

CoefficientImplementation::CoefficientImplementation(){}
CoefficientImplementation::~CoefficientImplementation(){}

bool CoefficientImplementation::equals(const Eigen::VectorXd& thisCoeff, 
                                       const Eigen::VectorXd& otherCoeff, 
                                       double tol) const {
  return (makeUniqueCopy(thisCoeff) - makeUniqueCopy(otherCoeff)).array().abs().maxCoeff() < tol;
}
 
Eigen::VectorXd CoefficientImplementation::makeUniqueCopy(const Eigen::VectorXd& thisCoeff) const {
  Eigen::VectorXd rval;
  makeUnique(thisCoeff, rval);
  return rval;
}

void CoefficientImplementation::print(const Eigen::VectorXd& thisCoeff, 
                                      const std::string& str) const {
  std::cout << str << " " << thisCoeff.transpose();
}


} // namespace curves
