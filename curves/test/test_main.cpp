/*
 * @file test_main.cpp
 * @date Aug 17, 2014
 * @author Paul Furgale
 */

#include <gtest/gtest.h>
#include <glog/logging.h>

/// Run all the tests that were declared with TEST()
int main(int argc, char **argv) {

  google::InitGoogleLogging(argv[0]);

  testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
