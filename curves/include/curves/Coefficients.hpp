#ifndef COEFFICIENTS_H_
#define COEFFICIENTS_H_

#include <string> // size_t

namespace curves {

class Coefficient;

/// \class Coefficients
/// \brief a pure virtual interface to a map of key to coefficient.
class Coefficients {
public:
  Coefficients();
  virtual ~Coefficients();

  /// \brief get the coefficient associated with the key.
  virtual const Coefficient& get(size_t key) const = 0;

private:
};

}  // namespace curves

#endif // COEFFICIENTS_H_
