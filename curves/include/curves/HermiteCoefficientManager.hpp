/*
 * @file HermiteCoefficientManager.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale, Abel Gawel, Renaud Dube
 */

#ifndef CT_HERMITE_COEFFICIENT_MANAGER_HPP
#define CT_HERMITE_COEFFICIENT_MANAGER_HPP

#include "CoefficientManagerBase.hpp"

namespace curves {

class HermiteCoefficientManager : public CoefficientManagerBase {
 public:
  HermiteCoefficientManager();
  virtual ~HermiteCoefficientManager();

  /// Print the value of the coefficient, for debugging and unit tests
  virtual void print(const std::string& str = "") const;

  /// \brief Get the coefficients that are active at a certain time.
  ///
  /// This method can fail if the time is out of bounds. If it
  /// Succeeds, the elements of the pair are guaranteed to be filled
  /// nonnull.
  ///
  /// @returns true if it was successful
  bool getCoefficientsAt(Time time, 
                         KeyCoefficientTime** outCoefficient0, KeyCoefficientTime** outCoefficient1) const;
};

} // namespace curves


#endif /* CT_HERMITE_COEFFICIENT_MANAGER_HPP */
