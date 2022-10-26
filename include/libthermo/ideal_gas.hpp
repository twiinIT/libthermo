// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_IDEAL_GAS_HPP
#define LIBTHERMO_IDEAL_GAS_HPP

#include "libthermo/thermo.hpp"
#include "libthermo/math_utils.hpp"

#include <cmath>


namespace thermo
{
    class IdealGas : public Thermo<IdealGas>
    {
    public:
        IdealGas(double r_, double cp_)
            : m_r(r_)
            , m_cp(cp_)
            , m_gamma(cp_ / (cp_ - r_)){};

        template <class P>
        IdealGas(const P& props)
            : m_r(props.r())
            , m_cp(props.cp())
            , m_gamma(props.cp() / (props.cp() - props.r())){};

        template <class T, IS_NOT_XTENSOR>
        auto gamma(const T& t) const;

        template <class T, IS_XTENSOR>
        auto gamma(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto cp(const T& = 0.) const;

        template <class T, IS_XTENSOR>
        auto cp(const T& = 0.) const;

        template <class T, IS_NOT_XTENSOR>
        auto phi(const T& t) const;

        template <class T, IS_XTENSOR>
        auto phi(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto pr(const T& t1, const T& t2, const T& eff_poly) const;

        template <class T, class E, IS_XTENSOR>
        auto pr(const T& t1, const T& t2, const E& eff_poly) const;

        template <class T, IS_NOT_XTENSOR>
        auto Tau(const T& p1, const T& p2, const T& eff_poly) const;

        template <class T, class E, IS_XTENSOR>
        auto Tau(const T& p1, const T& p2, const E& eff_poly) const;

        template <class T, IS_NOT_XTENSOR>
        auto eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T, IS_XTENSOR>
        auto eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T>
        auto h(const T& t) const;

        double r() const;

        double static_t(const double tt, const double mach, double = 1e-8, std::size_t = 1) const;

        template <class T>
        auto t_f_pr(
            const T& pr, const T& t1, const T& eff_poly, double = 1e-8, std::size_t = 1) const;

        template <class T>
        auto t_f_h(const T& h, double = 1e-8, std::size_t = 1) const;

        template <class T>
        auto t_f_phi(const T& h, double = 1e-8, std::size_t = 1) const;

        template <class T>
        auto mach_f_wqa(
            const T& pt, const T& tt, const T& wqa, double = 1e-8, std::size_t = 1) const;

    protected:
        double m_r, m_cp, m_gamma;
    };

    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas::gamma(const T&) const
    {
        return m_gamma;
    }

    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas::cp(const T&) const
    {
        return m_cp;
    }


    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas::phi(const T& t) const
    {
        return m_cp * std::log(t);
    }

    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas::pr(const T& t1, const T& t2, const T& eff_poly) const
    {
        return std::exp(std::log(t2 / t1) * eff_poly * m_cp / m_r);
    }

    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas::Tau(const T& p1, const T& p2, const T& eff_poly) const
    {
        return std::exp(std::log(p2 / p1) * m_r / (eff_poly * m_cp));
    }

    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas::eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return r() / m_cp * log(p2 / p1) / std::log(t2 / t1);
    }

    template <class T>
    auto IdealGas::h(const T& t) const
    {
        return m_cp * t;
    }

    inline double IdealGas::r() const
    {
        return m_r;
    }

    inline double IdealGas::static_t(const double tt, const double mach, double, std::size_t) const
    {
        return tt / (1 + 0.5 * (m_cp / (m_cp - m_r) - 1.) * std::pow(mach, 2.));
    }

    template <class T>
    auto IdealGas::t_f_pr(const T& pr, const T& t1, const T& eff_poly, double, std::size_t) const
    {
        return t1 * std::pow(pr, m_r / (m_cp * eff_poly));
    }

    template <class T>
    auto IdealGas::t_f_h(const T& h, double, std::size_t) const
    {
        return h / m_cp;
    }

    template <class T>
    auto IdealGas::t_f_phi(const T& phi, double, std::size_t) const
    {
        return std::exp(phi / m_cp);
    }

    template <class T>
    auto IdealGas::mach_f_wqa(const T& pt, const T& tt, const T& wqa, double, std::size_t) const
    {
        return 1.;
    }
}

// clang-format off
#ifdef LIBTHERMO_USE_XTENSOR
    #include "libthermo/detail/ideal_gas_xt_impl.hpp"
#endif

#endif
