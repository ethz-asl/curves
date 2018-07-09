cmake_minimum_required(VERSION 2.8.3)
project(curves)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

find_package(catkin QUIET)

if(catkin_FOUND)
  find_package(catkin REQUIRED COMPONENTS
  )
endif()

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

if(catkin_FOUND)
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
endif()

include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  ${kindr_INCLUDE_DIRS}
)
if(catkin_FOUND)
  include_directories(
    ${catkin_INCLUDE_DIRS}
  )
endif()

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
  ${Boost_LIBRARIES}
  glog
)
if(catkin_FOUND)
  target_link_libraries(${PROJECT_NAME}
    ${catkin_LIBRARIES}
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
endif()

if(catkin_FOUND)
  install(TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
    RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
  )

  install(DIRECTORY include/${PROJECT_NAME}/
    DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
  )
endif()

### Local Build ###

if(NOT catkin_FOUND)
  set(external_libs
    ${Boost_LIBRARIES}
    glog
  )

  # Register the local build in case on doesn't use "make install"
  export(PACKAGE curves)

  # Create variable with the library location
  set(LIB_NAME curves)
  set(curves_library $<TARGET_FILE:${LIB_NAME}>)

  # Create variable for the local build tree
  get_property(curves_include_dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
  # Configure config file for local build tree
  configure_file(curvesConfig.cmake.in
    "${PROJECT_BINARY_DIR}/curvesConfig.cmake" @ONLY
  )

### Install Build ###

  # Change the include location for the case of an install location
  set(curves_include_dirs
    ${CMAKE_INSTALL_PREFIX}/include
    ${EIGEN3_INCLUDE_DIR}
    ${kindr_INCLUDE_DIRS}
  )

  # We put the generated file for installation in a different repository (i.e., ./CMakeFiles/)
  configure_file(curvesConfig.cmake.in
    "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/curvesConfig.cmake" @ONLY
  )

  install(FILES
    "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/curvesConfig.cmake"
    DESTINATION share/curves/cmake COMPONENT dev
  )

  install(TARGETS ${PROJECT_NAME}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
  )

  install(DIRECTORY include/curves/
    DESTINATION include/curves
  )
endif()