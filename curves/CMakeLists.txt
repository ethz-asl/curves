cmake_minimum_required(VERSION 2.8.3)
project(curves)

find_package(catkin_simple 0.1.0 REQUIRED)

# Skip ALL_DEPS_REQUIRED so that we can run without gtsam.
catkin_simple()

option(USE_GTSAM "Should the curves library enable GTSAM-specific functionality?" 
	${gtsam_catkin_FOUND})

if(USE_GTSAM)
  message(STATUS "+++ Enabling GTSAM-specific functionality.")
  add_definitions(-DCURVES_USE_GTSAM)
endif()

cs_add_library(${PROJECT_NAME} 
  src/CoefficientImplementation.cpp
  src/VectorSpaceCoefficientImplementation.cpp
  src/CurveBase.cpp
  src/Coefficient.cpp
  src/KeyGenerator.cpp
  src/HermiteCoefficientManager.cpp
  src/EvaluatorBase.cpp
  src/Coefficients.cpp
  src/SlerpSE3Curve.cpp
  src/SE3CoefficientImplementation.cpp
  src/SE3Curve.cpp
)

find_package(Boost REQUIRED COMPONENTS system thread)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

catkin_add_gtest(${PROJECT_NAME}_tests 
  test/test_main.cpp
  test/test_HermiteCoefficientManager.cpp
  test/test_Coefficient.cpp
  test/test_Expressions.cpp
  test/test_Expressions_Laser.cpp
  test/test_Expressions_Slerp_SE3.cpp
  #test/test_SE3_Sinus_Circle.cpp
  test/test_MITb_dataset.cpp
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
)

message(glog_catkin_LIBRARIES: ${glog_catkin_LIBRARIES})
target_link_libraries(${PROJECT_NAME}_tests ${PROJECT_NAME})

add_doxygen(NOT_AUTOMATIC)

cs_install()
cs_export(CFG_EXTRAS curves-extras.cmake)

