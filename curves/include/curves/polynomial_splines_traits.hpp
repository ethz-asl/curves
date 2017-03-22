/*
 * polynomial_splines.hpp
 *
 *  Created on: Mar 7, 2017
 *      Author: dbellicoso
 */

#pragma once

// stl
#include <vector>

// eigen
#include <Eigen/Core>

// boost
#include <boost/math/special_functions/pow.hpp>


namespace curves {

struct SplineOpts {
  double tf = 0.0;
  double pos0 = 0.0;
  double posT = 1.0;
  double vel0 = 0.0;
  double velT = 0.0;
  double acc0 = 0.0;
  double accT = 0.0;
};

namespace spline_traits {

// Main struct template.
template<typename Core_, unsigned int SplineOrder_>
struct time_vector {
  const Eigen::Matrix<Core_, SplineOrder_+1, 1> tau(Core_ time);
  const Eigen::Matrix<Core_, SplineOrder_+1, 1> dtau(Core_ time);
  const Eigen::Matrix<Core_, SplineOrder_+1, 1> ddtau(Core_ time);

  const Eigen::Matrix<Core_, SplineOrder_+1, 1> tauZero();
  const Eigen::Matrix<Core_, SplineOrder_+1, 1> dtauZero();
  const Eigen::Matrix<Core_, SplineOrder_+1, 1> ddtauZero();

  bool compute(const SplineOpts& opts, std::vector<Core_>& coefficients);
};

// Specialization for third order splines.
template<>
struct time_vector<double, 3> {

  static constexpr unsigned int numCoefficients = 4;

  inline const std::vector<double> tau(double time) {
    return { boost::math::pow<3>(time),
             boost::math::pow<2>(time),
             time,
             1.0 };
  }

  inline const std::vector<double> dtau(double time) {
    return { 3.0*boost::math::pow<2>(time),
             2.0*time,
             1.0,
             0.0 };
  }

  inline const std::vector<double> ddtau(double time) {
    return { 6.0*time,
             2.0,
             0.0,
             0.0 };
  }

  inline const std::vector<double> tauZero()   { return { 0.0, 0.0, 0.0, 1.0 }; }
  inline const std::vector<double> dtauZero()  { return { 0.0, 0.0, 1.0, 0.0 }; }
  inline const std::vector<double> ddtauZero() { return { 0.0, 2.0, 0.0, 0.0 }; }

  bool compute(const SplineOpts& opts, std::vector<double>& coefficients) {
    const Eigen::Matrix<double, numCoefficients, 1> b = Eigen::Matrix<double, numCoefficients, 1>(opts.pos0, opts.vel0, opts.posT, opts.velT);
    Eigen::Matrix<double, numCoefficients, numCoefficients> A;
    A.row(0).data() = tauZero().data();
    A.row(1).data() = dtauZero().data();
    A.row(2).data() = tau(opts.tf).data();
    A.row(3).data() = dtau(opts.tf).data();
    const Eigen::VectorXd& coeffs = A.colPivHouseholderQr().solve(b);

    coefficients.data();

    for (int k=0; k<numCoefficients; k++) {
      coefficients[k] = coeffs(numCoefficients-1-k);
    }

    return true;
  }
};



// Specialization for fifth order splines.
template<>
struct time_vector<double, 5> {

  static constexpr unsigned int numCoefficients = 6;

  inline const std::vector<double> tau(double time) {
    return { boost::math::pow<5>(time),
             boost::math::pow<4>(time),
             boost::math::pow<3>(time),
             boost::math::pow<2>(time),
             time,
             1.0};
  }

  inline const std::vector<double> dtau(double time) {
    return { 5.0*boost::math::pow<4>(time),
             4.0*boost::math::pow<3>(time),
             3.0*boost::math::pow<2>(time),
             2.0*time,
             1.0,
             0.0};
  }

  inline const std::vector<double> ddtau(double time) {
    return { 20.0*boost::math::pow<3>(time),
             12.0*boost::math::pow<2>(time),
              6.0*time,
              2.0,
              0.0,
              0.0};
  }

  inline const std::vector<double> tauZero()   { return { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }; }
  inline const std::vector<double> dtauZero()  { return { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 }; }
  inline const std::vector<double> ddtauZero() { return { 0.0, 0.0, 0.0, 2.0, 0.0, 0.0 }; }

  bool compute(const SplineOpts& opts, std::vector<double>& coefficients) {
    const Eigen::Matrix<double, numCoefficients, 1> b = Eigen::Matrix<double, numCoefficients, 1>(opts.pos0, opts.vel0, opts.acc0, opts.posT, opts.velT, opts.accT);
    Eigen::Matrix<double, numCoefficients, numCoefficients> A;
    A.row(0) = Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>(tauZero());
    A.row(1) = Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>(dtauZero());
    A.row(2) = Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>(ddtauZero());
    A.row(3) = Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>(tau(opts.tf));
    A.row(4) = Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>(dtau(opts.tf));
    A.row(5) = Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>(ddtau(opts.tf));
    const Eigen::VectorXd& coeffs = A.colPivHouseholderQr().solve(b);

    coefficients.data();

    for (int k=0; k<numCoefficients; k++) {
      coefficients[k] = coeffs(5-k);
    }

    return true;
  }

};

}

}
