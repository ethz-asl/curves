#include <curves/Coefficient.hpp>
#include <curves/VectorSpaceCoefficientImplementation.hpp>
#include <glog/logging.h>

namespace curves {

Coefficient::Coefficient() : 
    impl_(new VectorSpaceCoefficientImplementation(0)) { }

Coefficient::Coefficient(size_t dim) : 
    impl_(new VectorSpaceCoefficientImplementation(dim)) { }

Coefficient::Coefficient(const Eigen::VectorXd& value) : 
    impl_(new VectorSpaceCoefficientImplementation(value.size())), 
    value_(value){ }


Coefficient::Coefficient(CoefficientImplementation::Ptr implementation) :
    impl_(implementation) {
  CHECK_NOTNULL(impl_.get());
  value_.resize(impl_->ambientDim());
}

Coefficient::Coefficient(CoefficientImplementation::Ptr implementation, const Eigen::VectorXd& value) :
    impl_(implementation), value_(value)
{
  CHECK_NOTNULL(impl_.get());
  CHECK_EQ(impl_->ambientDim(), value.size());
}

Coefficient::~Coefficient() {}

bool Coefficient::equals(const Coefficient& other, double tol) const {
  CHECK_EQ(value_.size(), other.value_.size());
  return impl_->equals(value_, other.value_, tol);
}
 
void Coefficient::print(const std::string& str) const {
  return impl_->print(value_, str);
}

size_t Coefficient::dim() const {
  return impl_->dim();
}

size_t Coefficient::ambientDim() const {
  return impl_->ambientDim();
}

Coefficient Coefficient::retract(const Eigen::VectorXd& delta) const {
  CHECK_EQ(delta.size(), impl_->dim());
  Coefficient rval(impl_);
  impl_->retract(value_, delta, rval.value_);
  return rval;
}

Eigen::VectorXd Coefficient::localCoordinates(const Coefficient& value) const {
  CHECK_EQ(value.value_.size(), impl_->ambientDim());
  Eigen::VectorXd rval;
  impl_->localCoordinates(value_, value.value_, rval);
  return rval;
}

const Eigen::VectorXd& Coefficient::getValue() const {
  return value_;
}

void Coefficient::setValue(const Eigen::VectorXd& value) {
  CHECK_EQ(value.size(), impl_->ambientDim());
  value_ = value;
}

double Coefficient::operator[](size_t i) const{
  CHECK_LT(i, value_.size());
  return value_[i];
}

} // namespace curves
