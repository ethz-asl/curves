/*
 * @file Support2CoefficientManager.hpp
 * @date Aug 17, 2014
 * @author Paul Furgale, Abel Gawel, Renaud Dube
 */

#ifndef CURVES_SUPPORT_2_COEFFICIENT_MANAGER_HPP
#define CURVES_SUPPORT_2_COEFFICIENT_MANAGER_HPP

#include "CoefficientManagerBase.hpp"

namespace curves {

class Support2CoefficientManager : public CoefficientManagerBase {
 public:
  Support2CoefficientManager();
  virtual ~Support2CoefficientManager();

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

#endif /* CURVES_SUPPORT_2_COEFFICIENT_MANAGER_HPP */
