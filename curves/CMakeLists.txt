cmake_minimum_required(VERSION 2.8.3)
project(curves)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

find_package(catkin REQUIRED COMPONENTS
  robot_utils
)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system thread)

# Attempt to find catkinized kindr
find_package(kindr QUIET)
if(NOT kindr_FOUND)
  # Attempt to find package-based kindr
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(kindr kindr REQUIRED)
endif()

# Add Doxygen documentation
add_subdirectory(doc/doxygen)

catkin_package(
  INCLUDE_DIRS
    include
  LIBRARIES
    ${PROJECT_NAME}
  CATKIN_DEPENDS
    robot_utils
  DEPENDS
    eigen
    boost
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${Kindr_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  src/KeyGenerator.cpp
  src/CubicHermiteSE3Curve.cpp
  src/CubicHermiteE3Curve.cpp
  src/SlerpSE3Curve.cpp
  src/SE3Curve.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${Boost_LIBRARIES}
)

add_dependencies(${PROJECT_NAME}
  ${catkin_EXPORTED_TARGETS}
)

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
target_link_libraries(${PROJECT_NAME}_tests
  ${PROJECT_NAME}
)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
