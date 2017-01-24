# Curves

A library for curves generation and estimation.

The source code is released under a [BSD 3-Clause license](ros_package_template/LICENSE).

**Authors: Renaud Dubé, Abel Gawel, Péter Fankhauser, Dario Bellicoso, Christian Gehring, Mike Bosse, Paul Furgale, Gabriel Agamennoni**

**Maintainer: Péter Fankhauser, pfankhauser@ethz.ch**

## Build status

[![Build Status](http://rsl-ci.ethz.ch/buildStatus/icon?job=curves)](http://rsl-ci.ethz.ch/job/curves/)

## Installation

### Installation from Packages

TODO

### Building from Source

#### Dependencies

- [Kindr](https://github.com/ethz-asl/kindr.git) (kinematics library)
- [Eigen](http://eigen.tuxfamily.org) (linear algebra library)

#### Building

To build from source, clone the latest version from this repository into your catkin workspace and compile the package using

	cd catkin_workspace/src
	git clone https://github.com/ethz-asl/curves.git
	cd ../
	catkin_make

### Unit Tests

Run the unit tests with

	catkin_make run_tests_curves run_tests_curves

## Bugs & Feature Requests

Please report bugs and request features using the [Issue Tracker](https://github.com/ethz-asl/curves/issues).
