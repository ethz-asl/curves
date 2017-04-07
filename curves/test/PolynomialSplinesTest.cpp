/*
 * PolynomialSplineContainerTest.cpp
 *
 *  Created on: Apr 2017
 *      Author: Dario Bellicoso
 *   Institute: ETH Zurich
 */

#include <gtest/gtest.h>

// curves
#include "curves/polynomial_splines.hpp"

// random number generation
#include <random>


// Construct a random number generator
std::random_device randomDevice;
std::default_random_engine randomEngine(randomDevice());
std::uniform_real_distribution<double> uniformDistribution(-10.0, 10.0);


TEST(PolynomialSplines, PolynomialSplinesCubic)
{
  curves::PolynomialSplineCubic spline;

  curves::SplineOptions opts(std::abs(uniformDistribution(randomEngine)),
                             uniformDistribution(randomEngine), uniformDistribution(randomEngine),
                             uniformDistribution(randomEngine), uniformDistribution(randomEngine),
                             uniformDistribution(randomEngine), uniformDistribution(randomEngine));

  spline.computeCoefficients(opts);

  EXPECT_NEAR(spline.getPositionAtTime(0.0), opts.pos0_, 1e-5);
  EXPECT_NEAR(spline.getPositionAtTime(opts.tf_), opts.posT_, 1e-5);

  EXPECT_NEAR(spline.getVelocityAtTime(0.0), opts.vel0_, 1e-5);
  EXPECT_NEAR(spline.getVelocityAtTime(opts.tf_), opts.velT_, 1e-5);
}

TEST(PolynomialSplines, PolynomialSplinesQuintic)
{
  curves::PolynomialSplineQuintic spline;

  curves::SplineOptions opts(std::abs(uniformDistribution(randomEngine)),
                             uniformDistribution(randomEngine), uniformDistribution(randomEngine),
                             uniformDistribution(randomEngine), uniformDistribution(randomEngine),
                             uniformDistribution(randomEngine), uniformDistribution(randomEngine));

  spline.computeCoefficients(opts);

  EXPECT_NEAR(spline.getPositionAtTime(0.0), opts.pos0_, 1e-5);
  EXPECT_NEAR(spline.getPositionAtTime(opts.tf_), opts.posT_, 1e-5);

  EXPECT_NEAR(spline.getVelocityAtTime(0.0), opts.vel0_, 1e-5);
  EXPECT_NEAR(spline.getVelocityAtTime(opts.tf_), opts.velT_, 1e-5);

  EXPECT_NEAR(spline.getAccelerationAtTime(0.0), opts.acc0_, 1e-5);
  EXPECT_NEAR(spline.getAccelerationAtTime(opts.tf_), opts.accT_, 1e-5);
}
