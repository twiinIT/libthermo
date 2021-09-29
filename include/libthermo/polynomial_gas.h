// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_POLY_GAS_H
#define LIBTHERMO_POLY_GAS_H

#include "libthermo/gas.h"

#include "xtensor/xcontainer.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xnoalias.hpp"

#include <string>

#include <xtensor/xfunction.hpp>



namespace libthermo
{
    using namespace xt;

    template<int D>
    struct poly
    {
        template<class E, class It>
        static inline auto polyval(E&& x, It coeff)
        {       
            return (*coeff++) + x * poly<D-1>::polyval(x, coeff);
        };
    };

    template<>
    struct poly<0>
    {
        template<class E, class It>
        static inline double polyval(E&& res, It coeff)
        {
            return 0;
        };
    };

    template<int I, class T, class It>
    static inline auto polyval(T&& t, It coeff)
    {
        return poly<I>::polyval(t, coeff);
    };

    class PolyGas : public Gas<PolyGas>
    {
    public:
        PolyGas(double r_, double cp_)
            : r(r_)
            , cp(cp_)
            , gamma(cp_ / (cp_ - r_))
        {};

        template<class T>
        auto Gamma(const T& t) const;

        template<class T>
        auto Cp(const T& t) const;

        template<class T>
        auto H(const T& t) const;

        template<class T, std::enable_if_t<xt::detail::is_container<T>::value, int> = 0>
        auto Phi(const T&t) const;

        template<class T, std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0>
        auto Phi(const T&t) const;

        template<class T, std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0>
        auto dPhi(const T& t1, const T& t2) const;

        template<class T, std::enable_if_t<xt::detail::is_container<T>::value, int> = 0>
        auto dPhi(const T& t1, const T& t2) const;

        template<class T>
        auto R() const;

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

        static const std::size_t coeff_size = 7;
        // https://www.cerfacs.fr/antares/_downloads/6d913fcaec57101ff421a2220b0769db/antares_doc.pdf p272
        const std::array<double, 8> cp_coeff = {2.35822e-20, -1.79613e-16, 4.70972e-13, -3.3094e-10, -6.27984e-07, 0.00123785, -0.521742, 1068.43};
        const std::array<double, 9> h_coeff = {2.35822e-20 / 8., -1.79613e-16 / 7., 4.70972e-13 / 6., -3.3094e-10 / 5., -6.27984e-07 / 4., 0.00123785 / 3., -0.521742 / 2., 1068.43, 0.};
        const std::array<double, 9> phi_coeff = {2.35822e-20 / 7., -1.79613e-16 / 6., 4.70972e-13 / 5., -3.3094e-10 / 4., -6.27984e-07 / 3., 0.00123785 / 2., -0.521742, 0., 1068.43};

    };

    template<class T>
    auto PolyGas::Gamma(const T& t) const
    {
        T tmp_cp = Cp(t);
        return tmp_cp / (tmp_cp - 287.05);
    }

    template<class T>
    auto PolyGas::Cp(const T& t) const
    {
        return polyval<7>(t, cp_coeff.rbegin());
    }

    template<class T>
    auto PolyGas::H(const T& t) const
    {
        return polyval<8>(t, h_coeff.rbegin());
    }

    template<class T, std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0>
    auto PolyGas::Phi(const T& t) const
    {
        return polyval<7>(t, ++phi_coeff.rbegin()) + phi_coeff[8] * std::log(t);
    }

    template<class T, std::enable_if_t<xt::detail::is_container<T>::value, int> = 0>
    auto PolyGas::Phi(const T& t) const
    {
        return polyval<7>(t, ++phi_coeff.rbegin()) + phi_coeff[8] * xt::log(t);
    }

    template<class T, std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0>
    auto PolyGas::dPhi(const T& t1, const T& t2) const
    {
        return (polyval<7>(t2, ++phi_coeff.rbegin()) - polyval<7>(t1, ++phi_coeff.rbegin())) + phi_coeff[8] * std::log(t2 / t1);
    }

    template<class T, std::enable_if_t<xt::detail::is_container<T>::value, int> = 0>
    auto PolyGas::dPhi(const T& t1, const T& t2) const
    {
        return (polyval<7>(t2, ++phi_coeff.rbegin()) - polyval<7>(t1, ++phi_coeff.rbegin())) + phi_coeff[8] * xt::log(t2 / t1);
    }

    template<class T>
    auto PolyGas::R() const
    {
        return r; 
    }

    template<class T>
    auto PolyGas::PR(const T& t1, const T& t2, const T& eff_poly) const
        -> std::enable_if_t<!xt::detail::is_container<T>::value, T>
    { 
        return std::exp(dPhi(t1, t2) * eff_poly / r);
    }

    template<class T, class E>
    auto PolyGas::PR(const T& t1, const T& t2, const E& eff_poly) const
        -> std::enable_if_t<xt::detail::is_container<T>::value, T>
    { 
        return xt::exp(dPhi(t1, t2) * eff_poly / r);
    }

    template<class T>
    auto PolyGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
        -> std::enable_if_t<!xt::detail::is_container<T>::value, T>
    {
        return R<T>() * log(p2 / p1) / dPhi(t1, t2);
    }

    template<class T>
    auto PolyGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
        -> std::enable_if_t<xt::detail::is_container<T>::value, T>
    {
        return R<double>() * xt::log(p2 / p1) / dPhi(t1, t2);
    }

    template<class T>
    T PolyGas::TFromH(const T& h) const
    { 
        return h / cp; 
    }
}
#endif
