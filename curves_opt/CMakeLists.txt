# Project configuration
cmake_minimum_required (VERSION 2.8)
project(curves_opt)

add_definitions(--std=c++11)

# Noisily default to Release build
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Find Catkin
find_package(catkin REQUIRED COMPONENTS
  curves
  numopt_quadprog
)

find_package(Eigen3)

# Attempt to find catkinized kindr
find_package(kindr QUIET)
if(NOT kindr_FOUND)
  # Attempt to find package-based kindr
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(kindr kindr REQUIRED)
endif()

# Declare this project as a catkin package
catkin_package(
  INCLUDE_DIRS
    include
    ${EIGEN3_INCLUDE_DIR}
  LIBRARIES
    ${PROJECT_NAME}
  CATKIN_DEPENDS
    curves
    numopt_quadprog
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${kindr_INCLUDE_DIRS}
)

add_library(
  ${PROJECT_NAME} 
  src/PolynomialSplineContainerOpt.cpp
)

target_link_libraries(
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
)

catkin_add_gtest(${PROJECT_NAME}_tests
  test/test_main.cpp
  test/PolynomialSplineContainerTest.cpp
  test/PolynomialSplineVectorSpaceCurveTest.cpp
  test/PolynomialSplineQuinticScalarCurveTest.cpp
  test/PolynomialSplinesTest.cpp
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/test
)

target_link_libraries(${PROJECT_NAME}_tests
  ${PROJECT_NAME}
  ${catkin_LIBRARIES}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
)

install(TARGETS ${PROJECT_NAME}
  ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
  RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
)
