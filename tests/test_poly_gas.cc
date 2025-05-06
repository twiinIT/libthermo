#include "libthermo/poly_gas.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "gtest/gtest.h"

#include <memory>
#include <iostream>


using namespace thermo;

TEST(PolyIdealGas, cp_tensor)
{
    auto props = PolyIdealGasProps<4>({ 0., 0., 0., 1004.4 }, 0., 287.05);
    auto dry_air = PolyIdealGas(props);

    EXPECT_EQ(dry_air.r(), 287.05);
    EXPECT_EQ(dry_air.h(1000.), 1004400.0);
    EXPECT_EQ(dry_air.cp(300.), 1004.4);
}
