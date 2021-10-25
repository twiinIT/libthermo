#include "libthermo/interp.h"

#include "xtensor/xtensor.hpp"

#include "gtest/gtest.h"


using namespace libthermo;


TEST(Interp1d, affine_scalar) 
{
  LinearFunc interpolator = { 1., 0.};

  EXPECT_DOUBLE_EQ(interpolator(1.), 1.);
  EXPECT_DOUBLE_EQ(interpolator(10.), 10.);
}


TEST(Interp1d, constant_scalar) 
{
  LinearFunc interpolator = { 0., 5.};

  EXPECT_DOUBLE_EQ(interpolator(1.), 5.);
  EXPECT_DOUBLE_EQ(interpolator(10.), 5.);
}

TEST(Interp1d, affine_tensor) 
{
  LinearFunc lin = { 2., 0.};
  const xt::xtensor<double, 1> t = { 0., 1., 10. };
  EXPECT_EQ(lin(t), 2 * t);
}


TEST(Interp1d, constant_tensor) 
{
  LinearFunc lin = { 0., 5.};
  const xt::xtensor<double, 1> t = { 0., 1., 10. };
  EXPECT_EQ(lin(t), xt::full_like(t, 5.));
}
