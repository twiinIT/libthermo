// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_IDEAL_GAS_H
#define LIBTHERMO_IDEAL_GAS_H

#ifdef LIBTHERMO_USE_XTENSOR
#include "xtensor/xcontainer.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xfunction.hpp"

#define IS_NOT_XTENSOR \
    std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0
#define IS_XTENSOR \
    std::enable_if_t<xt::detail::is_container<T>::value, int> = 0
#else
#include <type_traits>

namespace
{
    template<class T>
    struct is_tensor : std::false_type {};

    template<class T>
    struct is_not_tensor : std::true_type {};
}
#define IS_NOT_XTENSOR \
    std::enable_if_t<is_not_tensor<T>::value, int> = 0
#define IS_XTENSOR \
    std::enable_if_t<is_tensor<T>::value, int> = 0
#endif

#include "libthermo/gas.h"

#include <string>


namespace libthermo
{
    class IdealGas : public Gas<IdealGas>
    {
    public:
        IdealGas(double r_, double cp_)
            : r(r_)
            , cp(cp_)
            , gamma(cp_ / (cp_ - r_))
        {};

        template<class T, IS_NOT_XTENSOR>
        auto Gamma(const T& t) const;

        template<class T, IS_XTENSOR>
        auto Gamma(const T& t) const;

        template<class T, IS_NOT_XTENSOR>
        auto Cp(const T& = 0.) const;

        template<class T, IS_XTENSOR>
        auto Cp(const T& = 0.) const;

        template<class T, IS_NOT_XTENSOR>
        auto Phi(const T&t) const;

        template<class T, IS_XTENSOR>
        auto Phi(const T&t) const;

        template<class T, IS_NOT_XTENSOR>
        auto PR(const T& t1, const T& t2, const T& eff_poly) const;

        template<class T, class E, IS_XTENSOR>
        auto PR(const T& t1, const T& t2, const E& eff_poly) const;

        template<class T, IS_NOT_XTENSOR>
        auto EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template<class T, IS_XTENSOR>
        auto EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template<class T>
        auto H(const T& t) const;

        template<class T>
        auto R() const;

        template<class T>
        auto TFromH(const T& h) const;

    protected:
        double r, cp, gamma;
    };

    template<class T, IS_NOT_XTENSOR>
    auto IdealGas::Gamma(const T&) const
    {
        return gamma;
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template<class T, IS_XTENSOR>
    auto IdealGas::Gamma(const T& t) const
    {
        //return xt::full_like(t, gamma);
        //return xt::ones<double>(t.shape()) * gamma;
        return xt::broadcast(gamma, t.shape());
    }
#endif

    template<class T, IS_NOT_XTENSOR>
    auto IdealGas::Cp(const T&) const
    {
        return cp;
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template<class T, IS_XTENSOR>
    auto IdealGas::Cp(const T& t) const
    {
        return xt::full_like(t, cp);
    }
#endif

    template<class T, IS_NOT_XTENSOR>
    auto IdealGas::Phi(const T& t) const
    {
        return cp * std::log(t);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template<class T, IS_XTENSOR>
    auto IdealGas::Phi(const T& t) const
    {
        return cp * xt::log(t);
    }
#endif

    template<class T, IS_NOT_XTENSOR>
    auto IdealGas::PR(const T& t1, const T& t2, const T& eff_poly) const
    { 
        return std::exp(std::log(t2 / t1) * eff_poly * cp / r);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template<class T, class E, IS_XTENSOR>
    auto IdealGas::PR(const T& t1, const T& t2, const E& eff_poly) const
    { 
        return xt::exp(xt::log(t2 / t1) * eff_poly * cp / r);
    }
#endif

    template<class T, IS_NOT_XTENSOR>
    auto IdealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return R<T>() * cp * log(p2 / p1) / std::log(t2 / t1);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template<class T, IS_XTENSOR>
    auto IdealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return R<double>() * cp * xt::log(p2 / p1) / xt::log(t2 / t1);
    }
#endif

    template<class T>
    auto IdealGas::H(const T& t) const
    {
        return cp * t;
    }

    template<class T>
    auto IdealGas::R() const
    {
        return r; 
    }

    template<class T>
    auto IdealGas::TFromH(const T& h) const
    { 
        return h / cp; 
    }
}
#endif
