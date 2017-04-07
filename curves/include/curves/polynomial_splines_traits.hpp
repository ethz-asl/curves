/*
 * polynomial_splines.hpp
 *
 *  Created on: Mar 7, 2017
 *      Author: dbellicoso
 */

#pragma once

// stl
#include <vector>
#include <iostream>

// eigen
#include <Eigen/Core>
#include <Eigen/QR>

// boost
#include <boost/math/special_functions/pow.hpp>

namespace curves {

struct SplineOptions {

  constexpr SplineOptions(double tf, double pos0, double posT, double vel0,
                          double velT, double acc0, double accT)
      : tf_(tf),
        pos0_(pos0),
        posT_(posT),
        vel0_(vel0),
        velT_(velT),
        acc0_(acc0),
        accT_(accT)
  {

  }

  double tf_;
  double pos0_;
  double posT_;
  double vel0_;
  double velT_;
  double acc0_;
  double accT_;
};

namespace spline_traits {

// Main struct template.
template<typename Core_, int SplineOrder_>
struct spline_rep {

  static constexpr unsigned int splineOrder = SplineOrder_;
  static constexpr unsigned int numCoefficients = SplineOrder_+1;

  using TimeVectorType = std::array<Core_, numCoefficients>;
  using SplineCoefficients = std::array<double, numCoefficients>;

  static inline TimeVectorType tau(Core_ tk) noexcept;
  static inline TimeVectorType dtau(Core_ tk) noexcept;
  static inline TimeVectorType ddtau(Core_ tk) noexcept;

  static bool compute(const SplineOptions& opts, SplineCoefficients& coefficients);
};

// Specialization for third order splines.
template<>
struct spline_rep<double, 3> {

  static constexpr unsigned int splineOrder = 3;
  static constexpr unsigned int numCoefficients = splineOrder+1;

  using TimeVectorType = std::array<double, numCoefficients>;
  using SplineCoefficients = std::array<double, numCoefficients>;

  static inline TimeVectorType tau(double tk) noexcept {
    return { boost::math::pow<3>(tk), boost::math::pow<2>(tk), tk, 1.0 };
  }

  static inline TimeVectorType dtau(double tk) noexcept {
    return { 3.0*boost::math::pow<2>(tk), 2.0*tk, 1.0, 0.0 };
  }

  static inline TimeVectorType ddtau(double tk) noexcept {
    return { 6.0*tk, 2.0, 0.0, 0.0 };
  }

  static constexpr TimeVectorType   tauZero{{ 0.0, 0.0, 0.0, 1.0 }};
  static constexpr TimeVectorType  dtauZero{{ 0.0, 0.0, 1.0, 0.0 }};
  static constexpr TimeVectorType ddtauZero{{ 0.0, 2.0, 0.0, 0.0 }};

  static bool compute(const SplineOptions& opts, SplineCoefficients& coefficients) {

    Eigen::Matrix<double, numCoefficients, 1> b;
    b << opts.pos0_, opts.vel0_, opts.posT_, opts.velT_;

    Eigen::Matrix<double, numCoefficients, numCoefficients> A;
    A << Eigen::Map<const Eigen::Matrix<double, 1, numCoefficients>>(spline_rep<double, 3>::tauZero.data()),
         Eigen::Map<const Eigen::Matrix<double, 1, numCoefficients>>((dtau(0.0)).data()),
         Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>>((tau(opts.tf_)).data()),
         Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>>((dtau(opts.tf_)).data());

    const Eigen::VectorXd& coeffs = A.colPivHouseholderQr().solve(b);
    Eigen::Map<Eigen::VectorXd>( coefficients.data(), coeffs.rows(), coeffs.cols() ) = coeffs;

    return true;
  }
};


// Specialization for fifth order splines.
template<>
struct spline_rep<double, 5> {

  static constexpr unsigned int splineOrder = 5;
  static constexpr unsigned int numCoefficients = splineOrder+1;

  using TimeVectorType = std::array<double, numCoefficients>;
  using SplineCoefficients = std::array<double, numCoefficients>;

  static inline TimeVectorType tau(double tk) noexcept {
    return { boost::math::pow<5>(tk), boost::math::pow<4>(tk), boost::math::pow<3>(tk),
             boost::math::pow<2>(tk), tk, 1.0};
  }

  static inline TimeVectorType dtau(double tk) noexcept {
    return { 5.0*boost::math::pow<4>(tk), 4.0*boost::math::pow<3>(tk), 3.0*boost::math::pow<2>(tk),
             2.0*tk, 1.0, 0.0};
  }

  static inline TimeVectorType ddtau(double tk) noexcept {
    return { 20.0*boost::math::pow<3>(tk), 12.0*boost::math::pow<2>(tk), 6.0*tk,
              2.0, 0.0, 0.0};
  }

  static constexpr TimeVectorType   tauZero{{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 }};
  static constexpr TimeVectorType  dtauZero{{ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 }};
  static constexpr TimeVectorType ddtauZero{{ 0.0, 0.0, 0.0, 2.0, 0.0, 0.0 }};

  static bool compute(const SplineOptions& opts, SplineCoefficients& coefficients) {
    Eigen::Matrix<double, numCoefficients, 1> b;
    b << opts.pos0_, opts.vel0_, opts.acc0_, opts.posT_, opts.velT_, opts.accT_;

    Eigen::Matrix<double, numCoefficients, numCoefficients> A;
    A << Eigen::Map<const Eigen::Matrix<double, 1, numCoefficients>>(spline_rep<double, 5>::tauZero.data()),
         Eigen::Map<const Eigen::Matrix<double, 1, numCoefficients>>(spline_rep<double, 5>::dtauZero.data()),
         Eigen::Map<const Eigen::Matrix<double, 1, numCoefficients>>(spline_rep<double, 5>::ddtauZero.data()),
         Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>>((tau(opts.tf_)).data()),
         Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>>((dtau(opts.tf_)).data()),
         Eigen::Map<Eigen::Matrix<double, 1, numCoefficients>>((ddtau(opts.tf_)).data());

    const Eigen::VectorXd& coeffs = A.colPivHouseholderQr().solve(b);
    Eigen::Map<Eigen::VectorXd>( coefficients.data(), coeffs.rows(), coeffs.cols() ) = coeffs;

    return true;
  }

};

}

}
