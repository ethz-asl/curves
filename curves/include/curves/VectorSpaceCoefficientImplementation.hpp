#ifndef CT_VECTOR_SPACE_COEFFICIENT_IMPLEMENTATION_HPP
#define CT_VECTOR_SPACE_COEFFICIENT_IMPLEMENTATION_HPP

#include "CoefficientImplementation.hpp"

namespace curves {

class VectorSpaceCoefficientImplementation : public CoefficientImplementation
{
 public:
  /// \brief Initialize this with the dimension of
  ///        the vector space.
  VectorSpaceCoefficientImplementation(unsigned dimension);

  virtual ~VectorSpaceCoefficientImplementation();
  

  /// Compare this Coeficient with another for equality.
  virtual bool equals(const Eigen::VectorXd& thisCoeff, 
              const Eigen::VectorXd& otherCoeff, 
              double tol = 1e-9) const;

  /// For a given manifold, make the coefficient representation unique.
  /// This version modifies the argument in place
  virtual void makeUniqueInPlace(Eigen::VectorXd& thisCoeff) const;

  /// For a given manifold, make the coefficient representation unique.
  /// \param[in] thisCoeff The current value of the coefficient.
  /// \param[out] outUniqueCoeff The unique value of the coefficient.
  virtual void makeUnique(const Eigen::VectorXd& thisCoeff,
                  Eigen::VectorXd& outUniqueCoeff) const;
 
  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const Eigen::VectorXd& thisCoeff, 
             const std::string& str = "") const;

  /// Return the dimensionality of chosen chart for this coefficient.  This is
  /// the dimensionality of \c delta passed into retract() and of the vector
  /// returned by localCoordinates().
  /// @return The dimension of the chart
  virtual size_t dim() const;

  /// Return the dimensionality of the ambient space for this coefficient.  
  /// This is the dimensionality the actual coefficient vector.
  /// @return The dimension of the ambient space
  virtual size_t ambientDim() const;

  /// Increment the value, by mapping from the vector delta in the
  /// chart of the coefficient back to the manifold to produce a new,
  /// incremented value.
  /// @param[in] thisCoeff The current value of this coefficient
  /// @param[in] delta     The delta vector in the coefficient's chart, by
  ///                      which to increment this coefficient.
  /// @param[out] outIncrementedCoeff The output variable. 
  ///                      The incremented coefficient projected back
  ///                      to the manifold
  virtual void retract(const Eigen::VectorXd& thisCoeff, 
               const Eigen::VectorXd& delta, 
               Eigen::VectorXd& outIncrementedCoeff) const;

  /// Compute the coordinates in the chart assigned to this coefficient that
  /// retract() would map to \c value.
  /// @param value The value whose coordinates should be determined in this
  /// coefficient's chosen chart of the value on which this function is called.
  virtual void localCoordinates(const Eigen::VectorXd& thisCoeff, 
                        const Eigen::VectorXd& otherCoeff, 
                        Eigen::VectorXd& outLocalCoordinates) const;

 private:
  size_t dimension_;
};


} // namespace curves

#endif /* CT_VECTOR_SPACE_COEFFICIENT_IMPLEMENTATION_HPP */
