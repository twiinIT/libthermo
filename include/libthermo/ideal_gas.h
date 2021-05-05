// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_IDEAL_GAS_H
#define LIBTHERMO_IDEAL_GAS_H

#include "libthermo/gas.h"

#include "xtensor/xcontainer.hpp"
#include "xtensor/xtensor.hpp"

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

        template<class T, std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0>
        auto Gamma(const T& t) const;

        template<class T, std::enable_if_t<xt::detail::is_container<T>::value, int> = 0>
        auto Gamma(const T& t) const;

        template<class T>
        std::enable_if_t<!xt::detail::is_container<T>::value, T> Cp(const T& = 0.) const;

        template<class T>
        std::enable_if_t<xt::detail::is_container<T>::value, T> Cp(const T& = 0.) const;

        template<class T>
        T H(const T& t) const;

        template<class T>
        std::enable_if_t<!xt::detail::is_container<T>::value, T> Phi(const T&t) const;

        template<class T>
        std::enable_if_t<xt::detail::is_container<T>::value, T> Phi(const T&t) const;

        template<class T>
        T R() const;

        template<class T>
        std::enable_if_t<!xt::detail::is_container<T>::value, T> PR(const T& t1, const T& t2, const T& eff_poly) const;

        template<class T, class E>
        std::enable_if_t<xt::detail::is_container<T>::value, T> PR(const T& t1, const T& t2, const E& eff_poly) const;

        template<class T>
        std::enable_if_t<!xt::detail::is_container<T>::value, T> EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template<class T>
        std::enable_if_t<xt::detail::is_container<T>::value, T> EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template<class T>
        T TFromH(const T& h) const;

    protected:
        double r, cp, gamma;
    };

    template<class T, std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0>
    auto IdealGas::Gamma(const T&) const
    {
        return gamma;
    }

    template<class T, std::enable_if_t<xt::detail::is_container<T>::value, int> = 0>
    auto IdealGas::Gamma(const T& t) const
    {
        //return xt::full_like(t, gamma);
        //return xt::ones<double>(t.shape()) * gamma;
        return xt::broadcast(gamma, t.shape());
    }

    template<class T>
    auto IdealGas::Cp(const T&) const
        -> std::enable_if_t<!xt::detail::is_container<T>::value, T>
    {
        return cp;
    }

    template<class T>
    auto IdealGas::Cp(const T& t) const
        -> std::enable_if_t<xt::detail::is_container<T>::value, T>
    {
        return xt::full_like(t, cp);
    }

    template<class T>
    T IdealGas::H(const T& t) const
    {
        return cp * t;
    }

    template<class T>
    auto IdealGas::Phi(const T& t) const
        -> std::enable_if_t<!xt::detail::is_container<T>::value, T>
    {
        return cp * std::log(t);
    }

    template<class T>
    auto IdealGas::Phi(const T& t) const
        -> std::enable_if_t<xt::detail::is_container<T>::value, T>
    {
        return cp * xt::log(t);
    }

    template<class T>
    T IdealGas::R() const
    {
        return r; 
    }

    template<class T>
    auto IdealGas::PR(const T& t1, const T& t2, const T& eff_poly) const
        -> std::enable_if_t<!xt::detail::is_container<T>::value, T>
    { 
        return std::exp(std::log(t2 / t1) * eff_poly * cp / r);
    }

    template<class T, class E>
    auto IdealGas::PR(const T& t1, const T& t2, const E& eff_poly) const
        -> std::enable_if_t<xt::detail::is_container<T>::value, T>
    { 
        return xt::exp(xt::log(t2 / t1) * eff_poly * cp / r);
    }

    template<class T>
    auto IdealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
        -> std::enable_if_t<!xt::detail::is_container<T>::value, T>
    {
        return R<T>() * log(p2 / p1) / (Phi(t2) - Phi(t1));
    }

    template<class T>
    auto IdealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
        -> std::enable_if_t<xt::detail::is_container<T>::value, T>
    {
        return R<double>() * xt::log(p2 / p1) / (Phi(t2) - Phi(t1));
    }

    template<class T>
    T IdealGas::TFromH(const T& h) const
    { 
        return h / cp; 
    }
}
#endif
