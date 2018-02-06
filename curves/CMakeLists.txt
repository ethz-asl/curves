cmake_minimum_required(VERSION 2.8.3)
project(curves)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

find_package(catkin REQUIRED COMPONENTS
)

find_package(Eigen3 REQUIRED)
find_package(Boost REQUIRED COMPONENTS system thread)

# Glog
find_package(PkgConfig REQUIRED)
pkg_check_modules(glog libglog REQUIRED)

# Attempt to find catkinized kindr
find_package(kindr QUIET)
if(NOT kindr_FOUND)
  # Attempt to find package-based kindr
  pkg_check_modules(kindr kindr REQUIRED)
endif()

# Add Doxygen documentation
add_subdirectory(doc/doxygen)

catkin_package(
  INCLUDE_DIRS
    include
    ${EIGEN3_INCLUDE_DIR}
  LIBRARIES
    ${PROJECT_NAME}
  CATKIN_DEPENDS
  DEPENDS
    Boost
    glog
    kindr
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${kindr_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME}
  src/KeyGenerator.cpp
  src/CubicHermiteSE3Curve.cpp
  src/CubicHermiteE3Curve.cpp
  src/SlerpSE3Curve.cpp
  src/SE3Curve.cpp
  src/polynomial_splines_traits.cpp
  src/helpers.cpp
#  src/SE2Curve.cpp
#  src/SlerpSE2Curve.cpp
#  src/DiscreteSE3Curve.cpp
#  src/SemiDiscreteSE3Curve.cpp
#  src/SE3CurveFactory.cpp
)

target_link_libraries(${PROJECT_NAME}
  ${catkin_LIBRARIES}
  ${Boost_LIBRARIES}
  glog
)

add_dependencies(${PROJECT_NAME}
  ${catkin_EXPORTED_TARGETS}
)

catkin_add_gtest(${PROJECT_NAME}_tests
  test/test_main.cpp
  test/CubicHermiteSE3CurveTest.cpp
  test/PolynomialSplineContainerTest.cpp
  test/PolynomialSplineVectorSpaceCurveTest.cpp
  test/PolynomialSplineQuinticScalarCurveTest.cpp
  test/PolynomialSplinesTest.cpp
#  test/test_LocalSupport2CoefficientManager.cpp
#  test/test_Hermite.cpp
#  test/test_MITb_dataset.cpp
#  test/test_Pose2_Expressions.cpp
#  test/test_MITb_dataset_SE2.cpp
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
)

target_link_libraries(${PROJECT_NAME}_tests
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
  glog
)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)
