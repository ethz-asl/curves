# Curves

A library for curves generation and estimation.

The source code is released under a [BSD 3-Clause license](ros_package_template/LICENSE).

**Authors: Renaud Dubé, Abel Gawel, Péter Fankhauser, Dario Bellicoso, Christian Gehring, Mike Bosse, Paul Furgale, Gabriel Agamennoni**

**Maintainer: Péter Fankhauser, pfankhauser@anybotics.com**

## Build status

[![Build Status](https://ci.leggedrobotics.com/buildStatus/icon?job=github_ethz-asl/curves/master)](https://ci.leggedrobotics.com/job/github_ethz-asl/job/curves/job/master/)

## Installation

### Installation from Packages

TODO

### Building from Source

#### Dependencies

- [Eigen](http://eigen.tuxfamily.org) (linear algebra library)
- [Kindr](https://github.com/ethz-asl/kindr.git) (kinematics library)
- [Glog](https://github.com/google/glog) (logging library)

			sudo apt-get install libgoogle-glog-dev

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
