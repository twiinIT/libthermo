// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_POLY_GAS_HPP
#define LIBTHERMO_POLY_GAS_HPP

#include "libthermo/thermo.hpp"
#include "libthermo/poly_gas.hpp"
#include "libthermo/exceptions.hpp"
#include "libthermo/math_utils.hpp"
#include "libthermo/detail/polyval.hpp"

#include <boost/math/tools/roots.hpp>

#include <cmath>
#include <vector>
#include <exception>


namespace thermo
{
    template <class T, int D>
    struct PolyGasProps
    {
        using arr_t = std::array<T, D>;
        using value_t = T;

        PolyGasProps(const arr_t& cp, T h0, T rs);

        std::array<T, D> cp_coeffs;
        std::array<T, D - 1> dcp_coeffs;
        std::array<T, D + 1> h_coeffs;
        std::array<T, D> phi_coeffs;
        T r, phi_log;
    };

    template <class T, int D>
    PolyGasProps<T, D>::PolyGasProps(const arr_t& cp, T h0, T rs)
        : cp_coeffs(cp)
        , r(rs)
    {
        for (auto i = 0; i < D - 1; ++i)
        {
            dcp_coeffs[i] = cp_coeffs[i] * (D - i - 1);
            h_coeffs[i] = cp_coeffs[i] / (D - i);
            phi_coeffs[i] = cp_coeffs[i] / (D - i - 1);
        }
        h_coeffs[D - 1] = cp_coeffs[D - 1];
        h_coeffs[D] = h0;

        phi_coeffs[D - 1] = 0.;
        phi_log = cp_coeffs[D - 1];
    }

    template <class T, int D>
    PolyGasProps<T, D> mix(const std::vector<typename PolyGasProps<T, D>::arr_t>& cps,
                           const std::vector<T>& h0s,
                           const std::vector<T>& rs,
                           const std::vector<T>& weights)
    {
        typename PolyGasProps<T, D>::arr_t mix_cp;
        T mix_h0 = 0., mix_r = 0.;
        mix_cp.fill(0.);

        auto s = weights.size();
        if (s)
        {
            if (cps.size() != s || h0s.size() != s || rs.size() != s)
                throw std::runtime_error("Incorrect size");

            T total_weight = 0.;
            for (auto w : weights)
                total_weight += w;

            for (auto d = 0; d < D; ++d)
            {
                for (auto i = 0; i < s; ++i)
                    mix_cp[d] += cps[i][d] * weights[i];

                mix_cp[d] /= total_weight;
            }

            for (auto i = 0; i < s; ++i)
            {
                mix_h0 += h0s[i] * weights[i];
                mix_r += rs[i] * weights[i];
            }
            mix_h0 /= total_weight;
            mix_r /= total_weight;
        }

        return PolyGasProps<T, D>(mix_cp, mix_h0, mix_r);
    }

    template <class P>
    class PolyGas : public Thermo<PolyGas<P>>
    {
    public:
        using value_t = typename P::value_t;

        PolyGas(const P& properties)
            : m_props(properties){};

        template <class T, IS_NOT_XTENSOR>
        auto gamma(const T& t) const;

        template <class T, IS_XTENSOR>
        auto gamma(const T& t) const;

        template <class T>
        auto cp(const T& t) const;

        template <class T>
        auto h(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto phi(const T& t) const;

        template <class T, IS_XTENSOR>
        auto phi(const T& t) const;

        template <class T, IS_NOT_XTENSOR>
        auto dphi(const T& t1, const T& t2) const;

        template <class T, IS_XTENSOR>
        auto dphi(const T& t1, const T& t2) const;

        value_t r() const;

        template <class T, IS_NOT_XTENSOR>
        auto pr(const T& t1, const T& t2, const T& eff_poly) const;

        template <class T, class E, IS_XTENSOR>
        auto pr(const T& t1, const T& t2, const E& eff_poly) const;

        template <class T, IS_NOT_XTENSOR>
        auto eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T, IS_XTENSOR>
        auto eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const;

        template <class T>
        auto static_t(const T& tt, const T& mach, double tol, std::size_t max_iter = 30) const;

        template <class T>
        auto t_f_pr(const T& pr, const T& t1, const T& eff_poly, double tol, std::size_t max_iter = 30) const;

        template <class T>
        auto t_f_h(const T& h, double tol, std::size_t max_iter = 30) const;

        template <class T>
        auto t_f_phi(const T& h, double tol, std::size_t max_iter = 30) const;

        template <class T>
        auto mach_f_wqa(const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter = 30) const;

        const P& properties() const
        {
            return m_props;
        }

    protected:
        P m_props;
    };

    template <class P>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<P>::gamma(const T& t) const
    {
        auto tmp_cp = cp(t);
        return tmp_cp / (tmp_cp - m_props.r);
    }

    template <class P>
    template <class T>
    auto PolyGas<P>::cp(const T& t) const
    {
        using namespace detail;

        return polyval(t, m_props.cp_coeffs);
    }

    template <class P>
    template <class T>
    auto PolyGas<P>::h(const T& t) const
    {
        using namespace detail;

        return polyval(t, m_props.h_coeffs);
    }

    template <class P>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<P>::phi(const T& t) const
    {
        using namespace detail;

        return polyval(t, m_props.phi_coeffs) + m_props.phi_log * std::log(t);
    }

    template <class P>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<P>::dphi(const T& t1, const T& t2) const
    {
        using namespace detail;

        return (polyval(t2, m_props.phi_coeffs) - polyval(t1, m_props.phi_coeffs))
               + m_props.phi_log * std::log(t2 / t1);
    }

    template <class P>
    inline typename P::value_t PolyGas<P>::r() const
    {
        return m_props.r;
    }

    template <class P>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<P>::pr(const T& t1, const T& t2, const T& eff_poly) const
    {
        return std::exp(dphi(t1, t2) * eff_poly / m_props.r);
    }

    template <class P>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<P>::eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return r() * log(p2 / p1) / dphi(t1, t2);
    }

    template <class P>
    template <class T>
    inline auto PolyGas<P>::static_t(const T& tt, const T& mach, double tol, std::size_t max_iter) const
    {
        using namespace math;
        using namespace detail;

        double gam, dgam, ts, hs, v, cp_, dcp;
        bool converged = false;
        std::size_t iter;

        double ht = h(tt);
        double r_ = r();
        double x = 1.;
        double mach2 = pow(mach, 2.);

        // initialize using approximations
        gam = gamma(tt);
        ts = tt / (1. + (gam - 1.) / 2. * square(mach));

        for (auto i = 0; i < max_iter; ++i)
        {
            gam = gamma(ts);
            v = mach * std::sqrt(gam * r_ * ts);
            hs = h(ts);
            x = ht - hs - square(v) / 2.;

            if (std::abs(x / ht) < tol)
            {
                converged = true;
                iter = i;
                break;
            }

            dcp = polyval(ts, m_props.dcp_coeffs);
            cp_ = cp(ts);
            dgam = -r_ * dcp / pow(cp_ - r_, 2.);
            ts -= x / (-cp_ - 0.5 * mach2 * r_ * (gam + ts * dgam));
        }

        if (!converged)
            throw convergence_error();

        return ts;
    }

    template <class P>
    template <class T>
    auto PolyGas<P>::t_f_pr(const T& pr, const T& t1, const T& eff_poly, double tol, std::size_t max_iter) const
    {
        return t_f_phi(std::log(pr) * m_props.r / eff_poly + phi(t1), tol, max_iter);
    }

    template <class P>
    template <class T>
    auto PolyGas<P>::t_f_h(const T& h_in, double tol, std::size_t max_iter) const
    {
        using namespace detail;

        double t, cp_, h1, x;
        bool converged = false;

        t = 288.15;
        cp_ = cp(t);
        h1 = h(t);
        x = h_in - h1;

        for (auto i = 0; i < max_iter; ++i)
        {
            t += x / cp_;
            cp_ = cp(t);
            h1 = h(t);
            x = h_in - h1;

            if (std::abs(x / h_in) < tol)
            {
                converged = true;
                break;
            }
        }

        if (!converged)
            throw convergence_error();

        return t;
    }

    template <class P>
    template <class T>
    auto PolyGas<P>::t_f_phi(const T& phi_in, double tol, std::size_t max_iter) const
    {
        using namespace detail;

        double t, x1;
        bool converged = false;

        t = 288.15;
        x1 = phi_in - phi(t);

        for (auto i = 0; i < max_iter; ++i)
        {
            t += x1 * t / cp(t);
            x1 = phi_in - phi(t);

            if (std::abs(x1) < tol)
            {
                converged = true;
                break;
            }
        }

        if (!converged)
            throw convergence_error();

        return t;
    }

    template <class P>
    template <class T>
    auto PolyGas<P>::mach_f_wqa(const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter) const
    {
        using namespace detail;
        using namespace math;

        double ps, ts, ht, v, r_;
        bool converged = false;

        r_ = r();
        ht = h(tt);

        double ts_crit = static_t(tt, 1., tol);
        double ps_crit = pt * pr(tt, ts_crit, 1.);
        double wqa_crit = ps_crit / (r_ * ts_crit) * std::sqrt(gamma(ts_crit) * r_ * ts_crit);

        if (wqa < 0. || wqa > wqa_crit)
            throw domain_error();

        auto err_v = [&](double ts) -> double {
            ps = pt * pr(tt, ts, 1.);
            v = std::sqrt(2 * (ht - h(ts)));
            return ps / (r_ * ts) * v - wqa;
        };

        boost::uintmax_t niter = max_iter;
        auto res = boost::math::tools::toms748_solve(err_v,
                                                     ts_crit,
                                                     tt,
                                                     [&tol](const auto& a, const auto& b) -> bool {
                                                         using std::fabs;
                                                         return fabs(a - b) / (std::min)(fabs(a), fabs(b)) <= tol;
                                                     },
                                                     niter);
        ts = res.first;
        ps = pt * pr(tt, ts, 1.);
        v = std::sqrt(2 * (ht - h(ts)));

        return v / std::sqrt(gamma(ts) * r_ * ts);
    }
}

// clang-format off
#ifdef LIBTHERMO_USE_XTENSOR
    #include "libthermo/detail/poly_gas_xt_impl.hpp"
#endif

#endif
