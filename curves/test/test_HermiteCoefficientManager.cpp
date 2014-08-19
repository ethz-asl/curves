#include <gtest/gtest.h>
#include <curves/HermiteCoefficientManager.hpp>
#include <curves/VectorSpaceCoefficientImplementation.hpp>

class HermiteCoeffManagerTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    N = 50;
    for(size_t i = 0; i < N; ++i) {
      vectors.push_back(Eigen::VectorXd::Random(3));
      coefficients.push_back(curves::Coefficient(vectors[i]));
      // Make sure there are some negative times in there
      times.push_back(i * 1000 - 3250);
      keys1.push_back( manager1.insertCoefficient(times[i], coefficients[i]) );
    }
    manager2.insertCoefficients(times, coefficients, keys2);
  
    ASSERT_EQ(N, manager1.size());
    ASSERT_EQ(N, manager2.size());
  
    ASSERT_EQ(N, keys1.size());
    ASSERT_EQ(N, keys2.size());
  

    manager1.getTimes(times1);
    manager2.getTimes(times2);
  }

  // virtual void TearDown() {}
  
  size_t N;
  std::vector<Eigen::VectorXd> vectors;
  std::vector<curves::Coefficient> coefficients;
  std::vector<curves::Time> times;
  std::vector<curves::Time> times1;
  std::vector<curves::Time> times2;
  std::vector<curves::Key> keys1;
  std::vector<curves::Key> keys2;
  curves::HermiteCoefficientManager manager1;
  curves::HermiteCoefficientManager manager2;
  
};

TEST_F(HermiteCoeffManagerTest, testInsert) {
  
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


TEST_F(HermiteCoeffManagerTest, testTimes) {
  std::pair<curves::KeyCoefficientTime*, curves::KeyCoefficientTime*> bracket;
  bool success = false;
  curves::Time etime;
  
  etime = times[0] - 1;
  success = manager1.getCoefficientsAt(etime, bracket);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  etime = times[0] - 100;
  success = manager1.getCoefficientsAt(etime, bracket);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  etime = times[N-1];
  success = manager1.getCoefficientsAt(etime, bracket);
  ASSERT_TRUE(success) << "Eval at time " << etime;
  ASSERT_EQ(times[N-2],bracket.first->time) << "index " << N-2 << ", time: " << etime;
  ASSERT_EQ(times[N-1],bracket.second->time) << "index " << N-1 << ", time: " << etime;


  etime = times[N-1] + 1;
  success = manager1.getCoefficientsAt(etime, bracket);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  etime = times[N-1] + 100;
  success = manager1.getCoefficientsAt(etime, bracket);
  ASSERT_FALSE(success) << "Eval at time " << etime;

  for(size_t i = 1; i < times.size(); ++i) {

    etime = times[i-1];
    success = manager1.getCoefficientsAt(etime, bracket);
    ASSERT_TRUE(success) << "Eval at time " << etime;
    ASSERT_EQ(times[i-1],bracket.first->time) << "index " << i << ", time: " << etime;
    ASSERT_EQ(times[i],bracket.second->time) << "index " << i << ", time: " << etime;
    
    etime = (times[i-1] + times[i]) / 2;
    success = manager1.getCoefficientsAt(etime, bracket);
    ASSERT_TRUE(success) << "Eval at time " << etime;
    ASSERT_EQ(times[i-1],bracket.first->time) << "index " << i << ", time: " << etime;
    ASSERT_EQ(times[i],bracket.second->time) << "index " << i << ", time: " << etime;


  }
    
}


TEST_F(HermiteCoeffManagerTest, testGetCoefficientsInRange) {
  // \todo Abel and Renaud
}


TEST_F(HermiteCoeffManagerTest, testSetCoefficients) {
  // \todo Abel and Renaud
}

TEST_F(HermiteCoeffManagerTest, testRemoveCoefficients) {
  // \todo Abel and Renaud
}
