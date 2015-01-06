/*
 * @file CoefficientImplementation.cpp
 * @date Aug 17, 2014
 * @author Paul Furgale, Sean Anderson
 */

#ifndef CT_COEFFICIENT_IMPLEMENTATION_HPP
#define CT_COEFFICIENT_IMPLEMENTATION_HPP

#include <boost/shared_ptr.hpp>
#include <Eigen/Core>

namespace curves {

// \todo PTF Finish the documentation here
class CoefficientImplementation
{
 public:
  typedef boost::shared_ptr<CoefficientImplementation> Ptr;
  typedef boost::shared_ptr<const CoefficientImplementation> ConstPtr;

  CoefficientImplementation();
  virtual ~CoefficientImplementation();

  /// For a given manifold, make the coefficient representation unique.
  /// \param[in] thisCoeff The current value of the coefficient.
  /// \param[out] outUniqueCoeff The unique value of the coefficient.
  virtual void makeUnique(const Eigen::VectorXd& thisCoeff,
                          Eigen::VectorXd* outUniqueCoeff) const = 0;
 

  /// For a given manifold, make the coefficient representation unique.
  /// This version modifies the argument in place
  virtual void makeUniqueInPlace(Eigen::VectorXd* thisCoeff) const = 0;
  
 
  /// For a given manifold, make the coefficient representation unique.
  /// \param[in] thisCoeff The current value of the coefficient.
  /// \param[out] outUniqueCoeff The unique value of the coefficient.
  virtual Eigen::VectorXd makeUniqueCopy(const Eigen::VectorXd& thisCoeff) const;

  /// Project the coefficient back to the manifold. This is useful
  /// for manifolds that have constraints. For example, this will
  /// renormalize a unit-length quaternion coefficient that may have
  /// become unnormalized through repeated updates.
  //virtual Eigen::VectorXd projectToManifold(const Eigen::VectorXd& thisCoeff) const = 0;

  /// Project the coefficient back to the manifold. This is useful
  /// for manifolds that have constraints. For example, this will
  /// renormalize a unit-length quaternion coefficient that may have
  /// become unnormalized through repeated updates.
  //virtual Eigen::VectorXd& projectToManifoldInPlace(Eigen::VectorXd* thisCoeff) const = 0;


  /// Compare this Coefficient with another for equality.
  virtual bool equals(const Eigen::VectorXd& thisCoeff, 
                      const Eigen::VectorXd& otherCoeff, 
                      double tol = 1e-9) const;

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const Eigen::VectorXd& thisCoeff, 
                     const std::string& str = "") const;
 
  /// Return the dimensionality of chosen chart for this coefficient.  This is
  /// the dimensionality of \c delta passed into retract() and of the vector
  /// returned by localCoordinates().
  /// @return The dimension of the chart
  virtual size_t dim() const = 0;

  /// Return the dimensionality of the ambient space for this coefficient.  
  /// This is the dimensionality the actual coefficient vector.
  /// @return The dimension of the ambient space
  virtual size_t ambientDim() const = 0;

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
                       Eigen::VectorXd* outIncrementedCoeff) const = 0;

  /// Compute the coordinates in the chart assigned to this coefficient that
  /// retract() would map to \c value.
  /// @param value The value whose coordinates should be determined in this
  /// coefficient's chosen chart of the value on which this function is called.
  virtual void localCoordinates(const Eigen::VectorXd& thisCoeff, 
                                const Eigen::VectorXd& otherCoeff, 
                                Eigen::VectorXd* outLocalCoordinates) const = 0;
};

} // namespace curves

#endif /* CT_COEFFICIENT_IMPLEMENTATION_HPP */
