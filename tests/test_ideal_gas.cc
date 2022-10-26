#include "libthermo/ideal_gas.hpp"

#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "gtest/gtest.h"

#include <memory>


using namespace thermo;

TEST(IdealGas, cp_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.cp(288.15) << std::endl;

    xt::xtensor<double, 1> t = { 288.25, 400. };
    std::cout << gas.cp(t) << std::endl;

    xt::xtensor<double, 2> t2 = { { 288.25, 400. } };
    std::cout << gas.cp(t2) << std::endl;
}

TEST(IdealGas, h_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.h(288.15) << std::endl;

    xt::xtensor<double, 1> t = { 288.25, 400. };
    std::cout << gas.h(t) << std::endl;

    xt::xtensor<double, 2> t2 = { { 288.25, 400. } };
    std::cout << gas.h(t2) << std::endl;
}

TEST(IdealGas, phi_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.phi(288.15) << std::endl;

    xt::xtensor<double, 1> t = { 288.25, 400. };
    std::cout << gas.phi(t) << std::endl;

    xt::xtensor<double, 2> t2 = { { 288.25, 400. } };
    std::cout << gas.phi(t2) << std::endl;
}

TEST(IdealGas, pr_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.pr(288.15, 400., 0.8) << std::endl;
}

TEST(IdealGas, Constant)
{
    double constants[] = { 287, 290, 295 };
    for (auto r : constants)
    {
        IdealGas gas(r, 1004.685045);
        ASSERT_DOUBLE_EQ(gas.r(), r);
    }
}

TEST(IdealGas, gamma)
{
    double cps[] = { 1000., 1004, 1050 };
    double constants[] = { 287, 290, 295 };
    double temps[] = { 287., 500., 800., 1300., 1700. };
    for (auto cp : cps)
    {
        for (auto r : constants)
        {
            for (auto t : temps)
            {
                IdealGas gas(r, cp);
                ASSERT_DOUBLE_EQ(gas.gamma(t), cp / (cp - r));
            }
        }
    }
}

TEST(IdealGas, cp)
{
    double cps[] = { 1000., 1004, 1050 };
    double temps[] = { 287., 500., 800., 1300., 1700. };
    for (auto cp : cps)
    {
        for (auto t : temps)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.cp(t), cp);
        }
    }
}

TEST(IdealGas, h)
{
    double cps[] = { 1000., 1004, 1050 };
    double temps[] = { 287., 500., 800., 1300., 1700. };
    for (auto cp : cps)
    {
        for (auto t : temps)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.h(t), cp * t);
        }
    }
}

TEST(IdealGas, phi)
{
    double cps[] = { 1000., 1004, 1050 };
    double temps[] = { 287., 500., 800., 1300., 1700. };
    for (auto cp : cps)
    {
        for (auto t : temps)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.phi(t), cp * std::log(t));
        }
    }
}

TEST(IdealGas, t_f_h)
{
    double cps[] = { 1000., 1004, 1050 };
    double hs[] = { 287000., 500000., 800000., 1300000., 1700000. };
    for (auto cp : cps)
    {
        for (auto h : hs)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.t_f_h(h), h / cp);
        }
    }
}

TEST(IdealGas, pr)
{
    double taus[] = { 1., 1.2, 1.5, 3. };
    double eff_polys[] = { 0.2, 0.6, 1. };

    IdealGas gas(287., 1004.);

    for (auto tau : taus)
    {
        for (auto eff_poly : eff_polys)
        {
            ASSERT_DOUBLE_EQ(gas.pr(288.15, 288.15 * tau, eff_poly),
                             std::exp(std::log(tau) * eff_poly * gas.cp(288.15) / gas.r()));
        }
    }
}

TEST(IdealGas, eff_poly)
{
    double taus[] = { 1.15, 1.2, 1.5, 3. };
    double prs[] = { 1.2, 1.5, 1.7 };

    IdealGas gas(287., 1004.);
    double cp_ = gas.cp(288.15);

    for (auto tau : taus)
    {
        for (auto pr : prs)
        {
            ASSERT_DOUBLE_EQ(gas.eff_poly(101325., 288.15, 101325. * pr, 288.15 * tau),
                             gas.r() / cp_ * log(pr) / log(tau));
        }
    }
}
