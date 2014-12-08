/*
 * @file CoefficientManagerBase.hpp
 * @date Dec 3, 2014
 * @author Paul Furgale, Abel Gawel, Renaud Dube
 */

#ifndef CURVES_COEFFICIENT_MANAGER_BASE_HPP
#define CURVES_COEFFICIENT_MANAGER_BASE_HPP

#include <boost/unordered_map.hpp>
#include <vector>
#include <map>

#include "KeyCoefficientTime.hpp"

namespace curves {

class CoefficientManagerBase {
 public:
  CoefficientManagerBase();
  CoefficientManagerBase(const CoefficientManagerBase& rhs);
  CoefficientManagerBase& operator= (CoefficientManagerBase rhs);
  virtual ~CoefficientManagerBase();

  /// Compare this Coefficient with another for equality.
  virtual bool equals(const CoefficientManagerBase& other, double tol = 1e-9) const;

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// Get all of the keys in this manager. This method clears the
  /// list of keys before filling it.
  void getKeys(std::vector<Key>* outKeys) const;

  /// Get all of the keys in this manager. The list is not cleared
  /// before pushing it to the container.
  void appendKeys(std::vector<Key>* outKeys) const;

  /// Get a sorted list of coefficient times
  void getTimes(std::vector<Time>* outTimes) const;

  /// \brief insert a coefficient at a time and return
  ///        the key for the coefficient
  ///
  /// If a coefficient with this time already exists, it is overwritten
  Key insertCoefficient(Time time, const Coefficient& coefficient);

  /// \brief Insert coefficients. Optionally returns the keys for these coefficients.
  /// 
  /// If outKeys is not NULL, this function will not check if
  /// it is empty; new keys will be appended to this vector.
  void insertCoefficients(const std::vector<Time>& times,
                          const std::vector<Coefficient>& values,
                          std::vector<Key>* outKeys = NULL);

  /// \brief Remove the coefficient with this key.
  ///
  /// It is an error if the key does not exist.
  virtual void removeCoefficientWithKey(Key key);

  /// \brief Remove the coefficient at this time
  ///
  /// It is an error if there is no coefficient at this time.
  virtual void removeCoefficientAtTime(Time time);

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

  /// \brief Get the coefficients that are active within a range \f$[t_s,t_e) \f$.
  void getCoefficientsInRange(Time startTime, 
                              Time endTime, 
                              Coefficient::Map* outCoefficients) const;

  /// \brief Get all of the curve's coefficients.
  void getCoefficients(Coefficient::Map* outCoefficients) const;

  /// \brief Set coefficients.
  ///
  /// If any of these coefficients doen't exist, there is an error
  void setCoefficients(const Coefficient::Map& coefficients);

  /// \brief return the number of coefficients
  Key size() const;

  /// \brief clear the coefficients
  void clear();

  /// The first valid time for the curve.
  Time getMinTime() const;

  /// The one past the last valid time for the curve.
  Time getMaxTime() const;

  /// Check the internal consistency of the data structure
  /// If doExit is true, the function will call exit(0) at
  /// the end. This is useful for gtest death tests
  void checkInternalConsistency(bool doExit = false) const;

 protected:
  /// Key to coefficient mapping
  boost::unordered_map<Key, boost::shared_ptr<KeyCoefficientTime> > coefficients_;

  /// Time to coefficient mapping
  std::map<Time, boost::shared_ptr<KeyCoefficientTime> > timeToCoefficient_;

  /// \brief insert a coefficient that does not exist yet
  virtual void insertNewCoefficient(Key key, Time time, const Coefficient& coefficient);

  /// \brief check if a coefficient exists at time, return bool success and iterator
  bool hasCoefficientAtTime(Time time, std::map<Time, boost::shared_ptr<KeyCoefficientTime> >::iterator *it);

 private:
  /// Make new Coefficient container.
  /// Derived coefficient managers may override this
  /// function to use derivatives of KeyCoefficientTime
  virtual boost::shared_ptr<KeyCoefficientTime> makeContainer(Key key, Coefficient coefficient, Time time);
};

} // namespace curves

#endif /* CURVES_COEFFICIENT_MANAGER_BASE_HPP */
