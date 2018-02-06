/*
 * @file helpers.hpp
 * @date Oct 17, 2014
 * @author Sean Andersson, Peter Fankhauser
 */

#pragma once

#include <stdlib.h>
#include <vector>
#include <iostream>

#include <Eigen/Core>

#include "curves/Curve.hpp"

namespace curves {

template <typename T>
std::string toString(const T& t);

/// \brief Helper function to read CSV files into 'matrix' of strings
std::vector<std::vector<std::string> > loadCSV(std::string fileName);

/// \brief Helper function to write 'matrix' of strings into CSV file
void writeCSV(std::string fileName, const std::vector<std::vector<std::string> >& strMatrix);

/// \brief Helper function to read CSV files formatted in: time, vectorEntry0, vectorEntry1, ...
void loadTimeVectorCSV(std::string fileName, std::vector<curves::Time>* outTimes, std::vector<Eigen::VectorXd>* outValues);

/// \brief Helper function to write CSV file formatted in: time, vectorEntry0, vectorEntry1, ...
void writeTimeVectorCSV(std::string fileName, const std::vector<curves::Time>& times, const std::vector<Eigen::VectorXd>& values);

/// \brief Helper function to read CSV files formatted in: time0, time1, vectorEntry0, vectorEntry1, ...
void loadTimeTimeVectorCSV(std::string fileName, std::vector<curves::Time>* outTimes0, std::vector<curves::Time>* outTimes1, std::vector<Eigen::VectorXd>* outValues);

}
