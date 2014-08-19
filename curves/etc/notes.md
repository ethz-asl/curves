# Implementation Notes

* Should the *extend* methods and the *fit* methods return keys? Coefficients?

# Dealing with Values

* The big question, How to make this non-intrusive?
*  [This interface](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/base/DerivedValue.h?at=develop) seems tricky
, Specifically:
    * Supporting the boost::singleton_pool sounds like a tough job.
    * Figuring out a way to interface GTSAM without relying on GTSAM...
    * We have to make sure evaluators can be called functionally.
* Making coefficients be Eigen::VectorXd would be great.

## What are the requirements
* The Curve should
	* Have a unique key for each coefficient vector (Somebody has to explain how this works). This means it has a [Values](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/nonlinear/Values.h?at=develop) container for the whole spline
	* Return the vector of coefficients for each time
	* Return an evaluator for each time
	* **[Advanced]** Return evaluators that can move around in time (my gut tells me that if we get the Evaluator interface right, this can be added on without too much pain.
* The Evaluator should
	* implement something like `virtual Vector evaluate(const Values& x, boost::optional<std::vector<Matrix>&> H = boost::none) const;` to get the value and derivatives. This would be called from the Factor implementation. (still seems like this will be difficult).

# Dealing with Factors and Evaluators

Next: Look at the between factor to think about the implementation of the Evaluator.

From `NoiseModelFactor`

```c++
  /**
   * Error function *without* the NoiseModel, \f$ z-h(x) \f$.
   * Override this method to finish implementing an N-way factor.
   * If any of the optional Matrix reference arguments are specified, it should compute
   * both the function evaluation and its derivative(s) in X1 (and/or X2, X3...).
   */
  virtual Vector unwhitenedError(const Values& x, boost::optional< std::vector< Matrix > &> H = boost::none) const = 0;
```
So, Values is an key/Value map. The Factor decides the ordering that the Jacobians are returned.
# Curve
* `Curve::grow()` Needs to return a list of coefficients that are new
* `Curve::shrink()` needs to return a list of coefficients removed
# Open Questions
* Where do the keys live?
* How do we make factors that change the keys they are hooked up to dynamically (per iteration, for example)?
# Links
* [Documentation](https://research.cc.gatech.edu/borg/sites/edu.borg/files/downloads/gtsam.pdf)
## Factors
* [Factor](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/inference/Factor.h?at=develop)
* [GaussianFactor](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/linear/GaussianFactor.h?at=develop)
* [NonlinearFactor](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/nonlinear/NonlinearFactor.h?at=develop)
* [NoiseModelFactor](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/nonlinear/NonlinearFactor.h?at=develop)
## Values
* [Values](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/nonlinear/Values.h?at=develop)
* [Value](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/base/Value.h?at=develop) -- This has good documentation for writing a class derived from Value without using DerivedValue<T>. However, this documentation looks incomplete.
* [Derived Value](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/base/DerivedValue.h?at=develop)
* [LieVector](https://bitbucket.org/gtborg/gtsam/src/bcab483574f2bd636dcd8171cf1b62fdfe15b6d0/gtsam/base/LieVector.h?at=develop)
# Scratch
What about:

Coefficient
|
-- MyCoefficientType
|
-- MyGtsamCoefficientType : public MyCoefficientType, public DerivedValue<GtsamCoefficient>

* Coefficient implements equals, dim, retract, local coordinates, print...
* MyGtsamCoefficientType also implements these but forwards the implementation to the parent.

Is it better if we make it rely on GTSAM directly?
-- getCoefficients< CoefficientType* >
