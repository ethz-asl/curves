/*
 * @file test_LocalSupport2CoefficientManager.cpp
 * @date Aug 17, 2014
 * @author Paul Furgale
 */

#include <gtest/gtest.h>
#include <curves/LocalSupport2CoefficientManager.hpp>

using namespace curves;

class LocalSupport2CoefficientManagerTest : public ::testing::Test {
 protected:

  typedef Eigen::Matrix<double,3,1> Coefficient;
  typedef LocalSupport2CoefficientManager<Coefficient>::CoefficientIter CoefficientIter;

  virtual void SetUp() {
    N = 50;
    for(size_t i = 0; i < N; ++i) {
      coefficients.push_back(Coefficient::Random(3));
      // Make sure there are some negative times in there
      times.push_back(i * 1000 - 3250);
      keys1.push_back( manager1.insertCoefficient(times[i], coefficients[i]) );
    }
    manager2.insertCoefficients(times, coefficients, &keys2);

    ASSERT_EQ(N, manager1.size());
    ASSERT_EQ(N, manager2.size());

    ASSERT_EQ(N, keys1.size());
    ASSERT_EQ(N, keys2.size());

    manager1.getTimes(&times1);
    manager2.getTimes(&times2);
  }

  // virtual void TearDown() {}

  size_t N;
  std::vector<Coefficient> coefficients;
  std::vector<curves::Time> times;
  std::vector<curves::Time> times1;
  std::vector<curves::Time> times2;
  std::vector<curves::Key> keys1;
  std::vector<curves::Key> keys2;
  LocalSupport2CoefficientManager<Coefficient> manager1;
  LocalSupport2CoefficientManager<Coefficient> manager2;

};

TEST_F(LocalSupport2CoefficientManagerTest, testInsert) {

  ASSERT_EQ(N, manager1.size());
  ASSERT_EQ(N, manager2.size());

  ASSERT_EQ(N, keys1.size());
  ASSERT_EQ(N, keys2.size());

  ASSERT_EQ(N, times1.size());
  ASSERT_EQ(N, times2.size());

  for(size_t i = 0; i < N; ++i) {
    ASSERT_EQ(times1[i], times[i]);
    ASSERT_EQ(times2[i], times[i]);
  }

  ASSERT_EXIT(manager1.checkInternalConsistency(true), ::testing::ExitedWithCode(0), "^");
  ASSERT_EXIT(manager2.checkInternalConsistency(true), ::testing::ExitedWithCode(0), "^");
}


TEST_F(LocalSupport2CoefficientManagerTest, testTimes) {
  CoefficientIter bracket0;
  CoefficientIter bracket1;
  bool success = false;
  curves::Time etime;

  etime = times[0] - 1;
  success = manager1.getCoefficientsAt(etime, &bracket0, &bracket1);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  etime = times[0] - 100;
  success = manager1.getCoefficientsAt(etime, &bracket0, &bracket1);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  etime = times[N-1];
  success = manager1.getCoefficientsAt(etime, &bracket0, &bracket1);
  ASSERT_TRUE(success) << "Eval at time " << etime;
  ASSERT_EQ(times[N-2],bracket0->first) << "index " << N-2 << ", time: " << etime;
  ASSERT_EQ(times[N-1],bracket1->first) << "index " << N-1 << ", time: " << etime;

  etime = times[N-1] + 1;
  success = manager1.getCoefficientsAt(etime, &bracket0, &bracket1);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  etime = times[N-1] + 100;
  success = manager1.getCoefficientsAt(etime, &bracket0, &bracket1);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  for(size_t i = 1; i < times.size(); ++i) {

    etime = times[i-1];
    success = manager1.getCoefficientsAt(etime, &bracket0, &bracket1);
    ASSERT_TRUE(success) << "Eval at time " << etime;
    ASSERT_EQ(times[i-1],bracket0->first) << "index " << i << ", time: " << etime;
    ASSERT_EQ(times[i],bracket1->first) << "index " << i << ", time: " << etime;

    etime = (times[i-1] + times[i]) / 2;
    success = manager1.getCoefficientsAt(etime, &bracket0, &bracket1);
    ASSERT_TRUE(success) << "Eval at time " << etime;
    ASSERT_EQ(times[i-1],bracket0->first) << "index " << i << ", time: " << etime;
    ASSERT_EQ(times[i],bracket1->first) << "index " << i << ", time: " << etime;


  }

}

TEST_F(LocalSupport2CoefficientManagerTest, testGetCoefficientsInRange) {
  // \todo Abel and Renaud
}

TEST_F(LocalSupport2CoefficientManagerTest, testUpdateCoefficients) {

  for (size_t i = 0; i < keys1.size(); ++i) {
    manager1.updateCoefficientByKey(keys1[i], Coefficient::Zero());
    ASSERT_EQ(manager1.getCoefficientByKey(keys1[i]), Coefficient::Zero());
  }

  typedef boost::unordered_map<curves::Key, Coefficient> CoefficientMap;
  CoefficientMap allCoeffs;
  for (size_t i = 0; i < keys2.size(); ++i) {
    std::pair<curves::Key, Coefficient> pair = std::make_pair(keys2[i], Coefficient::Zero());
    allCoeffs.insert(pair);
  }

  manager2.updateCoefficients(allCoeffs);
  for (size_t i = 0; i < keys2.size(); ++i) {
    ASSERT_EQ(manager2.getCoefficientByKey(keys2[i]), Coefficient::Zero());
  }

}

TEST_F(LocalSupport2CoefficientManagerTest, testRemoveCoefficients) {
  // \todo Abel and Renaud
}
