cmake_minimum_required(VERSION 2.8.3)
project(curves)

find_package(catkin_simple REQUIRED)
catkin_simple(ALL_DEPS_REQUIRED)

add_definitions(-std=c++0x)

cs_add_library(${PROJECT_NAME} 
  src/CoefficientImplementation.cpp
  src/VectorSpaceCoefficientImplementation.cpp
  src/Curve.cpp
  src/Coefficient.cpp
  src/KeyGenerator.cpp
  src/HermiteCoefficientManager.cpp
  src/VectorSpaceCurve.cpp
  src/LinearInterpolationVectorSpaceCurve.cpp
  src/Evaluator.cpp
  src/VectorSpaceEvaluator.cpp
)

# find_package(Boost REQUIRED COMPONENTS system serialization filesystem)
# target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})

catkin_add_gtest(${PROJECT_NAME}_tests 
  test/test_main.cpp
  test/test_HermiteCoefficientManager.cpp
  test/test_Coefficient.cpp
)
target_link_libraries(${PROJECT_NAME}_tests ${PROJECT_NAME})

add_doxygen(REQUIRED)

cs_install()
cs_export()

