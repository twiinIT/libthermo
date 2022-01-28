// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_IDEAL_GAS_H
#define LIBTHERMO_IDEAL_GAS_H

#include "libthermo/gas.h"
#include "libthermo/math_utils.h"

#include <string>


namespace thermo
{
    class IdealGas : public Gas<IdealGas>
    {
    public:
        IdealGas(double r_, double cp_)
            : Gas<IdealGas>("IdealGas")
            , m_r(r_)
            , m_cp(cp_)
            , m_gamma(cp_ / (cp_ - r_)){};

        template <class T, IS_NOT_XTENSOR>
        auto Gamma(const T& t) const;

        template <class T, IS_XTENSOR>
        auto Gamma(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto Cp(const T& = 0.) const;

        template <class T, IS_XTENSOR>
        auto Cp(const T& = 0.) const;

        template <class T, IS_NOT_XTENSOR>
        auto Phi(const T& t) const;

        template <class T, IS_XTENSOR>
        auto Phi(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto PR(const T& t1, const T& t2, const T& eff_poly) const;

        template <class T, class E, IS_XTENSOR>
        auto PR(const T& t1, const T& t2, const E& eff_poly) const;

        template <class T, IS_NOT_XTENSOR>
        auto Tau(const T& p1, const T& p2, const T& eff_poly) const;

        template <class T, class E, IS_XTENSOR>
        auto Tau(const T& p1, const T& p2, const E& eff_poly) const;

        template <class T, IS_NOT_XTENSOR>
        auto EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T, IS_XTENSOR>
        auto EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T>
        auto H(const T& t) const;

        double R() const;

        double StaticT(const double tt, const double mach) const;

        template <class T>
        auto TFromPR(const T& pr, const T& t1, const T& eff_poly) const;

        template <class T>
        auto TFromH(const T& h) const;

        template <class T>
        auto TFromPhi(const T& h) const;

    protected:
        double m_r, m_cp, m_gamma;
    };

    template <class T, IS_NOT_XTENSOR>
    auto IdealGas::Gamma(const T&) const
    {
        return m_gamma;
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto IdealGas::Gamma(const T& t) const
    {
        // return xt::full_like(t, m_gamma);
        // return xt::ones<double>(t.shape()) * m_gamma;
        return xt::broadcast(m_gamma, t.shape());
    }
#endif

    template <class T, IS_NOT_XTENSOR>
    auto IdealGas::Cp(const T&) const
    {
        return m_cp;
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto IdealGas::Cp(const T& t) const
    {
        return xt::broadcast(m_cp, t.shape());
        // return xt::full_like(t, m_cp);
    }
#endif

    template <class T, IS_NOT_XTENSOR>
    auto IdealGas::Phi(const T& t) const
    {
        return m_cp * std::log(t);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto IdealGas::Phi(const T& t) const
    {
        return m_cp * xt::log(t);
    }
#endif

    template <class T, IS_NOT_XTENSOR>
    auto IdealGas::PR(const T& t1, const T& t2, const T& eff_poly) const
    {
        return std::exp(std::log(t2 / t1) * eff_poly * m_cp / m_r);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, class E, IS_XTENSOR>
    auto IdealGas::PR(const T& t1, const T& t2, const E& eff_poly) const
    {
        return xt::exp(xt::log(t2 / t1) * eff_poly * m_cp / m_r);
    }
#endif

    template <class T, IS_NOT_XTENSOR>
    auto IdealGas::Tau(const T& p1, const T& p2, const T& eff_poly) const
    {
        return std::exp(std::log(p2 / p1) * m_r / (eff_poly * m_cp));
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, class E, IS_XTENSOR>
    auto IdealGas::Tau(const T& p1, const T& p2, const E& eff_poly) const
    {
        return xt::exp(xt::log(p2 / p1) * m_r / (eff_poly * m_cp));
    }
#endif

    template <class T, IS_NOT_XTENSOR>
    auto IdealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return R() / m_cp * log(p2 / p1) / std::log(t2 / t1);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class T, IS_XTENSOR>
    auto IdealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return R() / m_cp * xt::log(p2 / p1) / xt::log(t2 / t1);
    }
#endif

    template <class T>
    auto IdealGas::H(const T& t) const
    {
        return m_cp * t;
    }

    inline double IdealGas::R() const
    {
        return m_r;
    }

    inline double IdealGas::StaticT(const double tt, const double mach) const
    {
        return tt / (1 + 0.5 * (m_cp / (m_cp - m_r) - 1.) * std::pow(mach, 2.));
    }

    template <class T>
    auto IdealGas::TFromPR(const T& pr, const T& t1, const T& eff_poly) const
    {
        return t1 * std::pow(pr, m_r / (m_cp * eff_poly));
    }

    template <class T>
    auto IdealGas::TFromH(const T& h) const
    {
        return h / m_cp;
    }

    template <class T>
    auto IdealGas::TFromPhi(const T& phi) const
    {
        return std::exp(phi / m_cp);
    }
}
#endif
