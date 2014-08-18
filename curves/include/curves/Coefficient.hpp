#ifndef CURVES_COEFFICIENT_HPP
#define CURVES_COEFFICIENT_HPP

#include <unordered_map>

#include <glog/logging.h>

#include "CoefficientImplementation.hpp"

namespace curves {

/// \class Coefficient
/// \brief A coefficient abstraction
///
/// \todo PTF shortly summarize the design
class Coefficient
{
 public:

  typedef std::unordered_map<size_t, Coefficient> Map;

  /// \brief The default constructor will make a zombie object that is not useful.
  Coefficient();

  /// \brief A convenience constructor to make a vector space coefficient.
  explicit Coefficient(size_t dim);

  /// \brief A convenience constructor to make a vector space coefficient.
  explicit Coefficient(const Eigen::VectorXd& v);

  /// \brief Initialize with just the implementation. The underlying vector
  ///        will be uninitialized.
  Coefficient(CoefficientImplementation::Ptr implementation);

  /// \brief The full constructor sets the implementation and the value 
  Coefficient(CoefficientImplementation::Ptr implementation, const Eigen::VectorXd& value);

  virtual ~Coefficient();

  /// Compare this Coeficient with another for equality.
  bool equals(const Coefficient& other, double tol = 1e-9) const;
 
  /// Print the value of the coefficient, for debugging and unit tests
  void print(const std::string& str = "") const;

  /// Return the dimensionality of chosen chart for this coefficient.  This is
  /// the dimensionality of \c delta passed into retract() and of the vector
  /// returned by localCoordinates().
  /// @return The dimension of the chart
  size_t dim() const;

  /// Return the dimensionality of the ambient space for this coefficient.  This is
  /// the dimensionality the actual coefficient vector.
  /// @return The dimension of the ambient space
  size_t ambientDim() const;

  /// Increment the value, by mapping from the vector delta in the
  /// chart of the coefficient back to the manifold to produce a new,
  /// incremented value.
  /// @param delta The delta vector in the coefficient's chart, by
  /// which to increment this coefficient.
  Coefficient retract(const Eigen::VectorXd& delta) const;

  /// Compute the coordinates in the chart assigned to this coefficient that
  /// retract() would map to \c value.
  /// @param value The value whose coordinates should be determined in this
  /// coefficient's chosen chart of the value on which this function is called.
  /// @return The coordinates of \c value in the chart.
  Eigen::VectorXd localCoordinates(const Coefficient& value) const;

  /// Get the value of the underlying coefficient
  const Eigen::VectorXd& getValue() const;

  /// Set the value of the underlying coefficient
  void setValue(const Eigen::VectorXd& value);

  /// Set the value of the underlying coefficient
  template<typename Derived>
  void setVector(const Eigen::MatrixBase<Derived>& value) {
    CHECK_EQ(value.rows(), ambientDim());
    CHECK_EQ(value.cols(), 1);
    value_ = value;
  }

  /// \brief get the value of the coefficient at element i
  double operator[](size_t i) const;
 private:
  Eigen::VectorXd value_;
  std::shared_ptr<CoefficientImplementation> impl_;

};


} // namespace curves

#endif /* CURVES_COEFFICIENT_HPP */
