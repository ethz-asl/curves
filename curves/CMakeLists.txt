cmake_minimum_required(VERSION 2.8.3)
project(curves)

find_package(catkin_simple 0.1.0 REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

set(CURVE_SRC
  include/curves/Coefficient.hpp
  include/curves/CoefficientImplementation.hpp
  include/curves/Coefficients.hpp
  include/curves/Curve.hpp
  include/curves/CurveBase.hpp
  include/curves/Evaluator.hpp
  include/curves/EvaluatorBase.hpp
  include/curves/GaussianProcessPriorFactorEvaluator.hpp
  include/curves/GaussianProcessVectorSpaceCurve.hpp
  include/curves/GaussianProcessVectorSpaceEvaluator.hpp
  include/curves/GaussianProcessVectorSpacePrior.hpp
  include/curves/HermiteCoefficientManager.hpp
  include/curves/KeyCoefficientTime.hpp
  include/curves/KeyGenerator.hpp
  include/curves/LinearInterpolationVectorSpaceCurve.hpp
  include/curves/LinearInterpolationVectorSpaceEvaluator.hpp
  include/curves/LinearSdeGaussianProcessPriorFactorEvaluator.hpp
  include/curves/LinearSdeGaussianProcessVectorSpacePrior.hpp
  include/curves/LtiCvGaussianProcessVectorSpaceCurve.hpp
  include/curves/LtiCvGaussianProcessVectorSpaceEvaluator.hpp
  include/curves/LtiCvGaussianProcessVectorSpacePrior.hpp
  include/curves/LtiSdeGaussianProcessVectorSpacePrior.hpp
  include/curves/ScalarCurveConfig.hpp
  include/curves/types.hpp
  include/curves/VectorSpaceCoefficientImplementation.hpp
  include/curves/VectorSpaceConfig.hpp
  include/curves/VectorSpaceCurve.hpp
  src/Coefficient.cpp
  src/Coefficients.cpp
  src/CoefficientImplementation.cpp
  src/CurveBase.cpp
  src/EvaluatorBase.cpp
  src/GaussianProcessPriorFactorEvaluator.cpp
  src/GaussianProcessVectorSpaceCurve.cpp
  src/GaussianProcessVectorSpaceEvaluator.cpp
  src/HermiteCoefficientManager.cpp
  src/KeyGenerator.cpp
  src/LinearInterpolationVectorSpaceCurve.cpp
  src/LinearInterpolationVectorSpaceEvaluator.cpp
  src/LinearSdeGaussianProcessPriorFactorEvaluator.cpp
  src/LinearSdeGaussianProcessVectorSpacePrior.cpp
  src/LtiCvGaussianProcessVectorSpaceCurve.cpp
  src/LtiCvGaussianProcessVectorSpaceEvaluator.cpp
  src/LtiCvGaussianProcessVectorSpacePrior.cpp
  src/LtiSdeGaussianProcessVectorSpacePrior.cpp
  src/VectorSpaceCoefficientImplementation.cpp
  src/VectorSpaceCurve.cpp
)

cs_add_library(${PROJECT_NAME} ${CURVE_SRC})

find_package(Boost REQUIRED COMPONENTS system thread)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

catkin_add_gtest(${PROJECT_NAME}_tests 
  test/test_main.cpp
  test/test_HermiteCoefficientManager.cpp
  test/test_Coefficient.cpp
)

message(glog_catkin_LIBRARIES: ${glog_catkin_LIBRARIES})
target_link_libraries(${PROJECT_NAME}_tests ${PROJECT_NAME})

add_doxygen(NOT_AUTOMATIC)

cs_install()
cs_export()

