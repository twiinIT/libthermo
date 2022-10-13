#include "libthermo/ideal_gas.h"

#include "xtensor/xtensor.hpp"
#include "xtensor/xio.hpp"

#include "gtest/gtest.h"

#include <memory>


using namespace thermo;

TEST(IdealGas, cp_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.Cp(288.15) << std::endl;

    xt::xtensor<double, 1> t = { 288.25, 400. };
    std::cout << gas.Cp(t) << std::endl;

    xt::xtensor<double, 2> t2 = { { 288.25, 400. } };
    std::cout << gas.Cp(t2) << std::endl;
}

TEST(IdealGas, h_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.H(288.15) << std::endl;

    xt::xtensor<double, 1> t = { 288.25, 400. };
    std::cout << gas.H(t) << std::endl;

    xt::xtensor<double, 2> t2 = { { 288.25, 400. } };
    std::cout << gas.H(t2) << std::endl;
}

TEST(IdealGas, phi_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.Phi(288.15) << std::endl;

    xt::xtensor<double, 1> t = { 288.25, 400. };
    std::cout << gas.Phi(t) << std::endl;

    xt::xtensor<double, 2> t2 = { { 288.25, 400. } };
    std::cout << gas.Phi(t2) << std::endl;
}

TEST(IdealGas, pr_tensor)
{
    IdealGas gas(287., 1004.685045);
    std::cout << gas.PR(288.15, 400., 0.8) << std::endl;
}

TEST(IdealGas, Constant)
{
    double constants[] = { 287, 290, 295 };
    for (auto r : constants)
    {
        IdealGas gas(r, 1004.685045);
        ASSERT_DOUBLE_EQ(gas.R(), r);
    }
}

TEST(IdealGas, Gamma)
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
                ASSERT_DOUBLE_EQ(gas.Gamma(t), cp / (cp - r));
            }
        }
    }
}

TEST(IdealGas, Cp)
{
    double cps[] = { 1000., 1004, 1050 };
    double temps[] = { 287., 500., 800., 1300., 1700. };
    for (auto cp : cps)
    {
        for (auto t : temps)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.Cp(t), cp);
        }
    }
}

TEST(IdealGas, H)
{
    double cps[] = { 1000., 1004, 1050 };
    double temps[] = { 287., 500., 800., 1300., 1700. };
    for (auto cp : cps)
    {
        for (auto t : temps)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.H(t), cp * t);
        }
    }
}

TEST(IdealGas, Phi)
{
    double cps[] = { 1000., 1004, 1050 };
    double temps[] = { 287., 500., 800., 1300., 1700. };
    for (auto cp : cps)
    {
        for (auto t : temps)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.Phi(t), cp * std::log(t));
        }
    }
}

TEST(IdealGas, TFromH)
{
    double cps[] = { 1000., 1004, 1050 };
    double hs[] = { 287000., 500000., 800000., 1300000., 1700000. };
    for (auto cp : cps)
    {
        for (auto h : hs)
        {
            IdealGas gas(287., cp);
            ASSERT_DOUBLE_EQ(gas.TFromH(h), h / cp);
        }
    }
}

TEST(IdealGas, PR)
{
    double taus[] = { 1., 1.2, 1.5, 3. };
    double eff_polys[] = { 0.2, 0.6, 1. };
    for (auto tau : taus)
    {
        for (auto eff_poly : eff_polys)
        {
            IdealGas gas(287., 1004.);
            ASSERT_DOUBLE_EQ(gas.PR(288.15, 288.15 * tau, eff_poly),
                             std::exp(std::log(tau) * eff_poly * gas.Cp(288.15) / gas.R()));
        }
    }
}

TEST(IdealGas, EffPoly)
{
    double taus[] = { 1.15, 1.2, 1.5, 3. };
    double prs[] = { 1.2, 1.5, 1.7 };
    for (auto tau : taus)
    {
        for (auto pr : prs)
        {
            IdealGas gas(287., 1004.);
            ASSERT_DOUBLE_EQ(gas.EffPoly(101325., 288.15, 101325. * pr, 288.15 * tau),
                             gas.R() * log(pr) / (gas.Phi(288.15 * tau) - gas.Phi(288.15)));
        }
    }
}
