// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_IDEAL_GAS_HPP
#define LIBTHERMO_IDEAL_GAS_HPP

#include "libthermo/thermo.hpp"
#include "libthermo/math_utils.hpp"

#include <boost/math/tools/roots.hpp>

#include <cmath>


namespace thermo
{
    template <class PropsType>
    class IdealGas : public Thermo<IdealGas<PropsType>>
    {
    public:
        IdealGas(const PropsType& r, const PropsType& cp)
            : m_r(r)
            , m_cp(cp)
            , m_gamma(cp / (cp - r)) {};

        template <class P>
        IdealGas(const P& props)
            : m_r(props.r())
            , m_cp(props.cp())
            , m_gamma(props.cp() / (props.cp() - props.r())) {};

        template <class T>
        auto gamma(const T& t) const;

        template <class T>
        auto cp(const T& = 0.) const;

        template <class T>
        auto phi(const T& t) const;

        template <class T>
        auto pr(const T& t1, const T& t2, const T& eff_poly) const;

        template <class T>
        auto Tau(const T& p1, const T& p2, const T& eff_poly) const;

        template <class T>
        auto eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T>
        auto h(const T& t) const;

        auto r() const;

        template <class T>
        auto static_t(const T& tt, const T& mach) const;
        template <class T>
        auto static_t(const T& tt, const T& mach, double, std::size_t = 1) const;

        template <class T>
        auto t_f_pr(const T& pr, const T& t1, const T& eff_poly) const;
        template <class T>
        auto t_f_pr(const T& pr, const T& t1, const T& eff_poly, double, std::size_t = 1) const;

        template <class T>
        auto t_f_h(const T& h) const;
        template <class T>
        auto t_f_h(const T& h, double, std::size_t = 1) const;

        template <class T>
        auto t_f_phi(const T& h) const;
        template <class T>
        auto t_f_phi(const T& h, double, std::size_t = 1) const;

        template <class T>
        auto mach_f_wqa(const T& pt, const T& tt, const T& wqa, double, std::size_t = 1) const;

    protected:
        PropsType m_r, m_cp, m_gamma;

        template <class S>
        friend bool operator==(const IdealGas<S>& ref, const IdealGas<S>& other);
    };

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::gamma(const T&) const
    {
        return m_gamma;
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::cp(const T&) const
    {
        return m_cp;
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::phi(const T& t) const
    {
        return m_cp * std::log(t);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::pr(const T& t1, const T& t2, const T& eff_poly) const
    {
        return std::exp(std::log(t2 / t1) * eff_poly * m_cp / m_r);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::Tau(const T& p1, const T& p2, const T& eff_poly) const
    {
        return std::exp(std::log(p2 / p1) * m_r / (eff_poly * m_cp));
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return r() / m_cp * log(p2 / p1) / std::log(t2 / t1);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::h(const T& t) const
    {
        return m_cp * t;
    }

    template <class PropsType>
    auto IdealGas<PropsType>::r() const
    {
        return m_r;
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::static_t(const T& tt, const T& mach) const
    {
        return tt / (1 + 0.5 * (m_cp / (m_cp - m_r) - 1.) * std::pow(mach, 2.));
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::static_t(const T& tt, const T& mach, double, std::size_t) const
    {
        return static_t(tt, mach);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::t_f_pr(const T& pr, const T& t1, const T& eff_poly) const
    {
        return t1 * std::pow(pr, m_r / (m_cp * eff_poly));
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::t_f_pr(const T& pr, const T& t1, const T& eff_poly, double, std::size_t) const
    {
        return t_f_pr(pr, t1, eff_poly);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::t_f_h(const T& h) const
    {
        return h / m_cp;
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::t_f_h(const T& h, double, std::size_t) const
    {
        return t_f_h(h);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::t_f_phi(const T& phi) const
    {
        return std::exp(phi / m_cp);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::t_f_phi(const T& phi, double, std::size_t) const
    {
        return t_f_phi(phi);
    }

    template <class PropsType>
    template <class T>
    auto IdealGas<PropsType>::mach_f_wqa(const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter) const
    {
        double a = wqa * sqrt(tt) / pt * sqrt(m_r / m_gamma);
        double coeff = 0.5 * (m_gamma + 1.) / (1. - m_gamma);
        double b = (m_gamma - 1.) / 2.;

        auto err_a = [&](double mach) -> double { return mach * pow(1. + b * mach * mach, coeff) - a; };

        boost::uintmax_t niter = max_iter;
        auto res = boost::math::tools::bracket_and_solve_root(
            err_a,
            0.5,
            1.2,
            true,
            [&tol](const auto& a, const auto& b) -> bool
            {
                using std::fabs;
                return fabs(a - b) / (std::min) (fabs(a), fabs(b)) <= tol;
            },
            niter);
        return res.first;
    }

    template <class PropsType>
    inline bool operator==(const IdealGas<PropsType>& ref, const IdealGas<PropsType>& other)
    {
        return ref.m_cp == other.m_cp && ref.m_r == other.m_r;
    }
}

#endif
