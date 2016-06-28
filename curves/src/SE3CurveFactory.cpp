/*
 * @file SE3CurveFactory.cpp
 * @date Feb 4th, 2016
 * @author Renaud Dubé, Abel Gawel
 */

#include <curves/SE3CurveFactory.hpp>

namespace curves {

std::shared_ptr<SE3Curve> SE3CurveFactory::create_curve(const std::string& curveType) {
  std::shared_ptr<SE3Curve> curve;

  std::cout << "Curve factory is creating a " << curveType << "." << std::endl;

  if(curveType == "slerp_curve") {
    curve = std::shared_ptr<SE3Curve>(new SlerpSE3Curve());
  } else if(curveType == "discrete_curve") {
    curve = std::shared_ptr<SE3Curve>(new DiscreteSE3Curve());
  } else if(curveType == "composition_curve") {
    curve = std::shared_ptr<SE3Curve>(new SE3CompositionCurve<SlerpSE3Curve, SlerpSE3Curve>());
  } else if(curveType == "discrete_composition_curve") {
    curve = std::shared_ptr<SE3Curve>(new SE3CompositionCurve<DiscreteSE3Curve, DiscreteSE3Curve>());
  } else if(curveType == "semi_discrete_composition_curve") {
    curve = std::shared_ptr<SE3Curve>(new SE3CompositionCurve<SemiDiscreteSE3Curve, SemiDiscreteSE3Curve>());
  } else {
    CHECK(false) << "This curve is not implemented.";
  }

  return curve;
}

}  // namespace curves
