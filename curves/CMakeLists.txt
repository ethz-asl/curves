cmake_minimum_required(VERSION 2.8.3)
project(curves)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

find_package(catkin_simple 0.1.0 REQUIRED)
find_package(Eigen3 REQUIRED)

# Skip ALL_DEPS_REQUIRED so that we can run without gtsam.
catkin_simple()

option(USE_GTSAM "Should the curves library enable GTSAM-specific functionality?" 
	${gtsam_catkin_FOUND})

if(USE_GTSAM)
  message(STATUS "+++ Enabling GTSAM-specific functionality.")
  add_definitions(-DCURVES_USE_GTSAM)
endif()

include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  ${Kindr_INCLUDE_DIRS}
)

# Attempt to find catkinized kindr
find_package(kindr QUIET)
if(NOT kindr_FOUND)
  # Attempt to find package-based kindr
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(kindr kindr REQUIRED)
endif()

cs_add_library(${PROJECT_NAME} 
  src/KeyGenerator.cpp
  src/CubicHermiteSE3Curve.cpp
  src/CubicHermiteE3Curve.cpp
  src/SlerpSE3Curve.cpp
  src/SE3Curve.cpp
)

find_package(Boost REQUIRED COMPONENTS system thread)
find_package(Eigen3 REQUIRED)
target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

catkin_add_gtest(${PROJECT_NAME}_tests 
  test/test_main.cpp
  test/CubicHermiteSE3CurveTest.cpp
  test/PolynomialSplineVectorSpaceCurveTest.cpp
  test/PolynomialSplineQuinticScalarCurveTest.cpp
#  test/test_LocalSupport2CoefficientManager.cpp
#  test/test_Expressions.cpp
#  test/test_Expressions_Slerp_SE3.cpp
#  test/test_Hermite.cpp
#  test/test_SE3_Sinus_Circle.cpp
#  test/test_SE3CompositionCurve.cpp
#  test/test_MITb_dataset.cpp
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
)

message(glog_catkin_LIBRARIES: ${glog_catkin_LIBRARIES})
target_link_libraries(${PROJECT_NAME}_tests ${PROJECT_NAME})

#add_doxygen(NOT_AUTOMATIC)

cs_install()
cs_export(CFG_EXTRAS curves-extras.cmake)

