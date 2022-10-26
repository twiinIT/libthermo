// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_DETAIL_REAL_GAS_XT_IMPL_HPP
#define LIBTHERMO_DETAIL_REAL_GAS_XT_IMPL_HPP

#include "libthermo/poly_gas.hpp"
#include "libthermo/type_traits.hpp"
#include "libthermo/detail/polyval.hpp"


namespace thermo
{
    template <class G>
    template <class T, IS_XTENSOR_>
    auto PolyGas<G>::gamma(const T& t) const
    {
        return cp(t) / (cp(t) - gas.r);
    }

    template <class G>
    template <class T, IS_XTENSOR_>
    auto PolyGas<G>::phi(const T& t) const
    {
        using namespace detail;
        return polyval(t, gas.phi_coeffs) + gas.phi_log * xt::log(t);
    }

    template <class G>
    template <class T, IS_XTENSOR_>
    auto PolyGas<G>::dphi(const T& t1, const T& t2) const
    {
        using namespace detail;
        return (polyval(t2, gas.phi_coeffs) - polyval(t1, gas.phi_coeffs))
               + gas.phi_log * xt::log(t2 / t1);
    }

    template <class G>
    template <class T, class E, IS_XTENSOR_>
    auto PolyGas<G>::pr(const T& t1, const T& t2, const E& eff_poly) const
    {
        return xt::exp(dphi(t1, t2) * eff_poly / gas.r);
    }

    template <class G>
    template <class T, IS_XTENSOR_>
    auto PolyGas<G>::eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return r() * xt::log(p2 / p1) / dphi(t1, t2);
    }
}
#endif
