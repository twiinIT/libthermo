// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_REAL_GAS_HPP
#define LIBTHERMO_REAL_GAS_HPP

#include "libthermo/thermo.hpp"

#ifdef LIBTHERMO_USE_XTENSOR
#include "xtensor/xexpression.hpp"
#include "xsimd/xsimd.hpp"
#define IS_NOT_XEXPRESSION std::enable_if_t<!xt::is_xexpression<E>::value, int> = 0
#define IS_XEXPRESSION std::enable_if_t<xt::is_xexpression<E>::value, int> = 0
#else
#include <cmath>
#include <type_traits>

namespace
{
    template <class T>
    struct is_xexpression : std::false_type
    {
    };

    template <class T>
    struct is_not_xexpression : std::true_type
    {
    };
}
#define IS_NOT_XEXPRESSION std::enable_if_t<is_not_tensor<E>::value, int> = 0
#define IS_XEXPRESSION std::enable_if_t<is_tensor<E>::value, int> = 0
#endif

#include <string>


namespace thermo
{
    template <int D>
    struct poly
    {
        template <class E, class It, IS_NOT_XEXPRESSION>
        static inline auto polyval(E&& x, It coeff)
        {
            auto c = coeff++;
            return std::fma(x, poly<D - 1>::polyval(x, coeff), *c);
        };

#ifdef LIBTHERMO_USE_XTENSOR
        template <class E, class It, IS_XEXPRESSION>
        static inline auto polyval(E&& x, It coeff)
        {
            auto c = coeff++;
            return xt::fma(x, poly<D - 1>::polyval(x, coeff), *c);
        };
#endif
    };

    template <>
    struct poly<0>
    {
        template <class E, class It>
        static inline double polyval(E&&, It coeff)
        {
            return *coeff;
        };
    };

    template <int I, class T, class It>
    static inline auto polyval(T&& x, It coeff)
    {
        return poly<I>::polyval(x, coeff);
    };

    class RealGas : public Thermo<RealGas>
    {
    public:
        RealGas(double r_)
            : Thermo<RealGas>("RealGas")
            , m_r(r_){};

        template <class T, IS_NOT_XTENSOR>
        auto Gamma(const T& t) const;

        template <class T, IS_XTENSOR>
        auto Gamma(const T& t) const;

        template <class T>
        auto Cp(const T& t) const;

        template <class T>
        auto H(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto Phi(const T& t) const;

        template <class T, IS_XTENSOR>
        auto Phi(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto dPhi(const T& t1, const T& t2) const;

        template <class T, IS_XTENSOR>
        auto dPhi(const T& t1, const T& t2) const;

        // template <class T>
        double R() const;

        template <class T, IS_NOT_XTENSOR>
        auto PR(const T& t1, const T& t2, const T& eff_poly) const;

        template <class T, class E, IS_XTENSOR>
        auto PR(const T& t1, const T& t2, const E& eff_poly) const;

        template <class T, IS_NOT_XTENSOR>
        auto EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T, IS_XTENSOR>
        auto EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T>
        auto StaticT(const T& tt, const T& mach) const;

        template <class T>
        auto TFromPR(const T& pr, const T& t1, const T& eff_poly) const;

        template <class T>
        auto TFromH(const T& h) const;

        template <class T>
        auto TFromPhi(const T& h) const;

    protected:
        double m_r;

        static const std::size_t coeff_size = 7;
        // https://www.cerfacs.fr/antares/_downloads/6d913fcaec57101ff421a2220b0769db/antares_doc.pdf
        // p272
        const std::array<double, 8> cp_coeff
            = { 2.35822e-20,  -1.79613e-16, 4.70972e-13, -3.3094e-10,
                -6.27984e-07, 0.00123785,   -0.521742,   1068.43 };
        const std::array<double, 9> h_coeff
            = { 2.35822e-20 / 8., -1.79613e-16 / 7., 4.70972e-13 / 6.,
                -3.3094e-10 / 5., -6.27984e-07 / 4., 0.00123785 / 3.,
                -0.521742 / 2.,   1068.43,           0. };
        const std::array<double, 9> phi_coeff = { 2.35822e-20 / 7.,
                                                  -1.79613e-16 / 6.,
                                                  4.70972e-13 / 5.,
                                                  -3.3094e-10 / 4.,
                                                  -6.27984e-07 / 3.,
                                                  0.00123785 / 2.,
                                                  -0.521742,
                                                  0.,
                                                  1068.43 };

        std::array<double, 5> reverse_coeffs
            = { 8.34098e-15, -2.36783e-10, 2.45527e-06, -0.010099, 18.2555 };
    };

    template <class T, IS_NOT_XTENSOR>
    auto RealGas::Gamma(const T& t) const
    {
        T tmp_cp = Cp(t);
        return tmp_cp / (tmp_cp - m_r);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto RealGas::Gamma(const T& t) const
    {
        return Cp(t) / (Cp(t) - m_r);
    }
#endif

    template <class T>
    auto RealGas::Cp(const T& t) const
    {
        return polyval<7>(t, cp_coeff.rbegin());
    }

    template <class T>
    auto RealGas::H(const T& t) const
    {
        return polyval<8>(t, h_coeff.rbegin());
    }

    template <class T, IS_NOT_XTENSOR>
    auto RealGas::Phi(const T& t) const
    {
        return polyval<7>(t, ++phi_coeff.rbegin()) + phi_coeff[8] * std::log(t);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto RealGas::Phi(const T& t) const
    {
        return polyval<7>(t, ++phi_coeff.rbegin()) + phi_coeff[8] * xt::log(t);
    }
#endif

    template <class T, IS_NOT_XTENSOR>
    auto RealGas::dPhi(const T& t1, const T& t2) const
    {
        return (polyval<7>(t2, ++phi_coeff.rbegin()) - polyval<7>(t1, ++phi_coeff.rbegin()))
               + phi_coeff[8] * std::log(t2 / t1);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto RealGas::dPhi(const T& t1, const T& t2) const
    {
        return (polyval<7>(t2, ++phi_coeff.rbegin()) - polyval<7>(t1, ++phi_coeff.rbegin()))
               + phi_coeff[8] * xt::log(t2 / t1);
    }
#endif

    // template <class T>
    inline double RealGas::R() const
    {
        return m_r;
    }

    template <class T, IS_NOT_XTENSOR>
    auto RealGas::PR(const T& t1, const T& t2, const T& eff_poly) const
    {
        return std::exp(dPhi(t1, t2) * eff_poly / m_r);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, class E, IS_XTENSOR>
    auto RealGas::PR(const T& t1, const T& t2, const E& eff_poly) const
    {
        return xt::exp(dPhi(t1, t2) * eff_poly / m_r);
    }
#endif

    template <class T, IS_NOT_XTENSOR>
    auto RealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return R() * log(p2 / p1) / dPhi(t1, t2);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto RealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return R() * xt::log(p2 / p1) / dPhi(t1, t2);
    }
#endif

    template <class T>
    inline auto RealGas::StaticT(const T& tt, const T& mach) const
    {
        return tt / (1 + 0.5 * (Cp(tt) / (Cp(tt) - m_r) - 1.) * std::pow(mach, 2.));
    }

    template <class T>
    auto RealGas::TFromPR(const T& pr, const T& t1, const T& eff_poly) const
    {
        return TFromPhi(std::log(pr) * m_r / eff_poly + Phi(t1));
    }

    template <class T>
    auto RealGas::TFromH(const T& h) const
    {
        double t, cp, h1, x;
        std::size_t counter = 0;
        std::array<double, 7> reverse_coeffs_h
            = { 4.94485e-36, -4.36014e-29, 1.46887e-22, -2.19316e-16,
                6.79214e-11, 0.0010072,    -10.4407 };

        t = polyval<6>(h, reverse_coeffs_h.rbegin());
        cp = Cp(t);
        h1 = H(t);
        x = h - h1;

        while (counter < 16 && std::abs(x) > 1e-10)
        {
            t += x / cp;
            cp = Cp(t);
            h1 = H(t);
            x = h - h1;
            ++counter;
        }

        return t;
    }

    template <class T>
    auto RealGas::TFromPhi(const T& phi) const
    {
        double t, cp, x, x1;
        std::size_t counter = 0;

        x = polyval<4>(phi, reverse_coeffs.rbegin());
        t = std::exp(x);

        x1 = phi - Phi(t);
        while (counter < 20 && std::abs(x1) > 1e-10)
        {
            t += x1 * t / Cp(t);
            x1 = phi - Phi(t);
            ++counter;
        }

        return t;
    }
}
#endif
