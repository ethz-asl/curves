#include <curves/LinearInterpolationVectorSpaceEvaluator.hpp>
#include <iostream>
using namespace std;

namespace curves {

void LinearInterpolationVectorSpaceEvaluator::getKeys(std::vector<Key> *outKeys) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

void LinearInterpolationVectorSpaceEvaluator::getCoefficients(std::vector<Coefficient> *outCoefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluate(const std::vector<Coefficient>& coefficients) const {
  cout << "You evaluated an object!" << endl;
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluateAndJacobian(const std::vector<Coefficient>& coefficients,
                                                                                                                std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the ambient space of the curve.
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateVector() const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives.
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the ambient space of the curve (functional form).
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateVector(const std::vector<Coefficient>& coefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives (functional form).
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateDerivative(unsigned derivativeOrder,
                                                                            const std::vector<Coefficient>& coefficients) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the ambient space of the curve (functional form).
LinearInterpolationVectorSpaceEvaluator::ValueType LinearInterpolationVectorSpaceEvaluator::evaluateVectorAndJacobian(const std::vector<Coefficient>& coefficients,
                                                                                                                      std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

/// Evaluate the curve derivatives (functional form).
Eigen::VectorXd LinearInterpolationVectorSpaceEvaluator::evaluateDerivativeAndJacobian(unsigned derivativeOrder,
                                                                                       const std::vector<Coefficient>& coefficients,
                                                                                       std::vector<Eigen::MatrixXd>* outJacobian) const {
  // \todo Abel and Renaud
  CHECK(false) << "To be implemented by Abel and Renaud :-)";
}

} // namespace curves
