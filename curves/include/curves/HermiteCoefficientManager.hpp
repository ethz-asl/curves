#ifndef CT_HERMITE_COEFFICIENT_MANAGER_HPP
#define CT_HERMITE_COEFFICIENT_MANAGER_HPP

#include <unordered_map>
#include <vector>
#include <map>

#include "KeyCoefficientTime.hpp"

namespace curves {

class HermiteCoefficientManager {
 public:
  HermiteCoefficientManager();
  virtual ~HermiteCoefficientManager();

  /// Compare this Coeficient with another for equality.
  bool equals(const HermiteCoefficientManager& other, double tol = 1e-9) const;
 
  /// Print the value of the coefficient, for debugging and unit tests
  void print(const std::string& str = "") const;

  /// Get all of the keys in this manager. This method clears the
  /// list of keys before filling it.
  void getKeys(std::vector<Key>& outKeys) const;

  /// Get all of the keys in this manager. The list is not cleared
  /// before pushing it to the container.
  void appendKeys(std::vector<Key>& outKeys) const;

  /// Get a sorted list of coefficient times
  void getTimes(std::vector<Time>& outTimes) const;
  
  /// \brief insert a coefficient at a time and return
  ///        the key for the coefficient
  ///
  /// If a coefficient with this time already exists, it is overwritten
  Key insertCoefficient(Time time, const Coefficient& coefficient);

  /// \brief Insert coefficients. Returns the keys for these coefficients.
  /// 
  /// This function will not check if outKeys is empty. New keys will
  /// be appended to this vector.
  void insertCoefficients(const std::vector<Time>& times,
                          const std::vector<Coefficient>& values,
                          std::vector<Key>& outKeys);
  
  /// \brief Remove the coefficient with this key.
  ///
  /// It is an error if the key does not exist.
  void removeCoefficientWithKey(Key key);

  /// \brief Remove the coefficient at this time
  ///
  /// It is an error if there is no coefficient at this time.
  void removeCoefficientAtTime(Time time);

  /// \brief return true if there is a coefficient at this time
  bool hasCoefficientAtTime(Time time) const;

  /// \brief return true if there is a coefficient with this key
  bool hasCoefficientWithKey(Key key) const;

  /// \brief set the coefficient associated with this key
  ///
  /// This function fails if there is no coefficient associated
  /// with this key.
  void setCoefficientByKey(Key key, const Coefficient& coefficient);

  /// \brief set the coefficient associated with this key
  ///
  /// This function fails if there is no coefficient associated
  /// with this key.
  void setCoefficientVectorByKey(Key key, const Eigen::VectorXd& vector);

  /// \brief get the coefficient associated with this key
  Coefficient getCoefficientByKey(Key key) const;
  
  /// \brief Get the coefficients that are active at a certain time.
  ///
  /// This method can fail if the time is out of bounds. If it
  /// Succeeds, the elements of the pair are guaranteed to be filled
  /// nonnull.
  ///
  /// @returns true if it was successful
  bool getCoefficientsAt(Time time, 
                         std::pair<KeyCoefficientTime*, KeyCoefficientTime*>& outCoefficients) const;

  
  /// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
  void getCoefficientsInRange(Time startTime, 
                              Time endTime, 
                              Coefficient::Map& outCoefficients) const;
  
  /// \brief Get all of the curve's coefficients.
  void getCoefficients(Coefficient::Map& outCoefficients) const;

  /// \brief Set coefficients.
  ///
  /// If any of these coefficients doen't exist, there is an error
  void setCoefficients(Coefficient::Map& coefficients);

  /// \brief return the number of coefficients
  Key size() const;
  
  /// \brief clear the coefficients
  void clear();

  /// The first valid time for the curve.
  Time getBackTime() const;
  
  /// The one past the last valid time for the curve.
  Time getFrontTime() const;

  /// Check the internal consistency of the data structure
  /// If doExit is true, the function will call exit(0) at
  /// the end. This is useful for gtest death tests
  void checkInternalConsistency(bool doExit = false) const;
  
 private:
  /// Key to coefficient mapping
  std::unordered_map<Key, KeyCoefficientTime> coefficients_;
  
  /// Time to coefficient mapping
  std::map<Time, KeyCoefficientTime*> timeToCoefficient_;
};

} // namespace curves


#endif /* CT_HERMITE_COEFFICIENT_MANAGER_HPP */
