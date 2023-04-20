// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_DETAIL_IDEAL_GAS_XT_IMPL_HPP
#define LIBTHERMO_DETAIL_IDEAL_GAS_XT_IMPL_HPP

#include "libthermo/ideal_gas.hpp"
#include "libthermo/type_traits.hpp"


namespace thermo
{
    template <class D>
    template <class T, IS_XTENSOR_>
    auto IdealGas<D>::gamma(const T& t) const
    {
        return xt::broadcast(m_gamma, t.shape());
    }

    template <class D>
    template <class T, IS_XTENSOR_>
    auto IdealGas<D>::cp(const T& t) const
    {
        return xt::broadcast(m_cp, t.shape());
    }

    template <class D>
    template <class T, IS_XTENSOR_>
    auto IdealGas<D>::phi(const T& t) const
    {
        return m_cp * xt::log(t);
    }

    template <class D>
    template <class T, class E, IS_XTENSOR_>
    auto IdealGas<D>::pr(const T& t1, const T& t2, const E& eff_poly) const
    {
        return xt::exp(xt::log(t2 / t1) * eff_poly * m_cp / m_r);
    }

    template <class D>
    template <class T, class E, IS_XTENSOR_>
    auto IdealGas<D>::Tau(const T& p1, const T& p2, const E& eff_poly) const
    {
        return xt::exp(xt::log(p2 / p1) * m_r / (eff_poly * m_cp));
    }

    template <class D>
    template <class T, IS_XTENSOR_>
    auto IdealGas<D>::eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return r() / m_cp * xt::log(p2 / p1) / xt::log(t2 / t1);
    }
}
#endif
