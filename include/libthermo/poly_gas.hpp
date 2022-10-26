// Copyright (c) 2021-2022, twiinIT
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

#include <vector>
#include <exception>


namespace thermo
{
    template <int D>
    struct PolyGasProps
    {
        using arr_t = std::array<double, D>;

        PolyGasProps(const arr_t& cp_coeffs_, double h0_, double r_);

        std::array<double, D> cp_coeffs;
        std::array<double, D + 1> h_coeffs;
        std::array<double, D> phi_coeffs;
        double r, phi_log;
    };

    template <int D>
    PolyGasProps<D>::PolyGasProps(const arr_t& cp_coeffs_, double h0_, double r_)
        : cp_coeffs(cp_coeffs_)
        , r(r_)
    {
        for (auto i = 0; i < D - 1; ++i)
        {
            h_coeffs[i] = cp_coeffs[i] / (D - i);
            phi_coeffs[i] = cp_coeffs[i] / (D - i - 1);
        }
        h_coeffs[D - 1] = cp_coeffs[D - 1];
        h_coeffs[D] = h0_;

        phi_coeffs[D - 1] = 0.;
        phi_log = cp_coeffs[D - 1];
    }

    template <int D>
    PolyGasProps<D> mix(const std::vector<typename PolyGasProps<D>::arr_t>& cps_,
                        const std::vector<double>& h0s_,
                        const std::vector<double>& rs_,
                        const std::vector<double>& weights_)
    {
        typename PolyGasProps<D>::arr_t mix_cp;
        double mix_h0 = 0., mix_r = 0.;
        mix_cp.fill(0.);

        auto s = weights_.size();
        if (s)
        {
            if (cps_.size() != s || h0s_.size() != s || rs_.size() != s)
                throw std::runtime_error("Incorrect size");

            double total_weight = 0.;
            for (auto w : weights_)
                total_weight += w;

            for (auto d = 0; d < D; ++d)
            {
                for (auto i = 0; i < s; ++i)
                    mix_cp[d] += cps_[i][d] * weights_[i];

                mix_cp[d] /= total_weight;
            }

            for (auto i = 0; i < s; ++i)
            {
                mix_h0 += h0s_[i] * weights_[i];
                mix_r += rs_[i] * weights_[i];
            }
            mix_h0 /= total_weight;
            mix_r /= total_weight;
        }

        return PolyGasProps<D>(mix_cp, mix_h0, mix_r);
    }

    template <class D>
    class PolyGas : public Thermo<PolyGas<D>>
    {
    public:
        PolyGas(const D& gas_)
            : gas(gas_){};

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

        double r() const;

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
        auto t_f_pr(const T& pr,
                    const T& t1,
                    const T& eff_poly,
                    double tol,
                    std::size_t max_iter = 30) const;

        template <class T>
        auto t_f_h(const T& h, double tol, std::size_t max_iter = 30) const;

        template <class T>
        auto t_f_phi(const T& h, double tol, std::size_t max_iter = 30) const;

        template <class T>
        auto mach_f_wqa(
            const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter = 30) const;

    protected:
        D gas;
    };

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<D>::gamma(const T& t) const
    {
        T tmp_cp = cp(t);
        return tmp_cp / (tmp_cp - gas.r);
    }

    template <class D>
    template <class T>
    auto PolyGas<D>::cp(const T& t) const
    {
        using namespace detail;

        return polyval(t, gas.cp_coeffs);
    }

    template <class D>
    template <class T>
    auto PolyGas<D>::h(const T& t) const
    {
        using namespace detail;

        return polyval(t, gas.h_coeffs);
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<D>::phi(const T& t) const
    {
        using namespace detail;

        return polyval(t, gas.phi_coeffs) + gas.phi_log * std::log(t);
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<D>::dphi(const T& t1, const T& t2) const
    {
        using namespace detail;

        return (polyval(t2, gas.phi_coeffs) - polyval(t1, gas.phi_coeffs))
               + gas.phi_log * std::log(t2 / t1);
    }

    template <class D>
    inline double PolyGas<D>::r() const
    {
        return gas.r;
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<D>::pr(const T& t1, const T& t2, const T& eff_poly) const
    {
        return std::exp(dphi(t1, t2) * eff_poly / gas.r);
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto PolyGas<D>::eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return r() * log(p2 / p1) / dphi(t1, t2);
    }

    template <class D>
    template <class T>
    inline auto PolyGas<D>::static_t(const T& tt,
                                     const T& mach,
                                     double tol,
                                     std::size_t max_iter) const
    {
        using namespace math;

        double gam, ts, hs, v, cp_;
        bool converged = false;

        double ht = h(tt);
        double r_ = r();
        double x = 1.;

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
                break;
            }

            cp_ = cp(ts);
            ts -= x / (-cp_ - v * mach * std::sqrt(gam * r_) / (2. * std::sqrt(ts)));
        }

        if (!converged)
            throw convergence_error();

        return ts;
    }

    template <class D>
    template <class T>
    auto PolyGas<D>::t_f_pr(
        const T& pr, const T& t1, const T& eff_poly, double tol, std::size_t max_iter) const
    {
        return t_f_phi(std::log(pr) * gas.r / eff_poly + phi(t1), tol, max_iter);
    }

    template <class D>
    template <class T>
    auto PolyGas<D>::t_f_h(const T& h_in, double tol, std::size_t max_iter) const
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

            if (std::abs(x) < tol)
            {
                converged = true;
                break;
            }
        }

        if (!converged)
            throw convergence_error();

        return t;
    }

    template <class D>
    template <class T>
    auto PolyGas<D>::t_f_phi(const T& phi_in, double tol, std::size_t max_iter) const
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

    template <class D>
    template <class T>
    auto PolyGas<D>::mach_f_wqa(
        const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter) const
    {
        using namespace detail;
        using namespace math;

        double ps, ts, hs, ht, v, r_, cp_, x, der, rho_v;
        bool converged = false;

        r_ = r();
        ht = h(tt);

        ts = static_t(tt, 0.5, tol);
        x = 1.;

        for (auto i = 0; i < max_iter; ++i)
        {
            hs = h(ts);
            v = std::sqrt(2 * (ht - hs));
            ps = pr(tt, ts, 1.) * pt;
            rho_v = ps / (ts * r_) * v;

            x = rho_v / wqa - 1.;
            if (std::abs(x) < tol)
            {
                converged = true;
                break;
            }

            cp_ = cp(ts);
            der = (rho_v / wqa) * (cp_ * (1.0 / (r_ * ts) - 1.0 / square(v)) - 1.0 / ts);
            ts -= x / der;
        }

        if (!converged)
            throw convergence_error();

        return v / std::sqrt(gamma(ts) * r_ * ts);
    }
}

// clang-format off
#ifdef LIBTHERMO_USE_XTENSOR
    #include "libthermo/detail/poly_gas_xt_impl.hpp"
#endif

#endif
