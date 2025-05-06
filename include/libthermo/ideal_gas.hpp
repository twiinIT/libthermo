// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_IDEAL_GAS_HPP
#define LIBTHERMO_IDEAL_GAS_HPP

#include "libthermo/thermo.hpp"
#include "libthermo/math_utils.hpp"
#include "libthermo/thread_pool.hpp"

#include <boost/math/tools/roots.hpp>

#include <xsimd/xsimd.hpp>

#include <cmath>


namespace thermo
{
    template <class D>
    class IdealGas : public Thermo<IdealGas<D>>
    {
    public:
        IdealGas(const D& r, const D& cp)
            : m_r(r)
            , m_cp(cp)
            , m_gamma(cp / (cp - r)) {};

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
        void eff_poly(T* res, T* p1, T* t1, T* p2, T* t2, std::size_t size) const;

        template <class T>
        auto h(const T& t) const;

        template <class T>
        void h(T* res, T* temperature, std::size_t size) const;

        double r() const;

        double static_t(const double tt, const double mach) const;
        double static_t(const double tt, const double mach, double, std::size_t = 1) const;

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
        auto mach_f_wqa(const T& pt, const T& tt, const T& wqa) const;
        template <class T>
        auto mach_f_wqa(const T& pt, const T& tt, const T& wqa, double, std::size_t = 1) const;

    protected:
        D m_r, m_cp, m_gamma;

        template <class S>
        friend bool operator==(const IdealGas<S>& ref, const IdealGas<S>& other);
    };

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas<D>::gamma(const T&) const
    {
        return m_gamma;
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas<D>::cp(const T&) const
    {
        return m_cp;
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas<D>::phi(const T& t) const
    {
        return m_cp * std::log(t);
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas<D>::pr(const T& t1, const T& t2, const T& eff_poly) const
    {
        return std::exp(std::log(t2 / t1) * eff_poly * m_cp / m_r);
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas<D>::Tau(const T& p1, const T& p2, const T& eff_poly) const
    {
        return std::exp(std::log(p2 / p1) * m_r / (eff_poly * m_cp));
    }

    template <class D>
    template <class T, IS_NOT_XTENSOR_>
    auto IdealGas<D>::eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return r() / m_cp * std::log(p2 / p1) / std::log(t2 / t1);
    }


    template <class D>
    template <class T>
    void IdealGas<D>::eff_poly(T* res, T* p1, T* t1, T* p2, T* t2, std::size_t size) const
    {
        using b_type = xsimd::batch<double, xsimd::avx2>;
        std::size_t inc = b_type::size;

        std::size_t vec_size = size - size % inc;
        b_type cp = m_cp;

        for (std::size_t i = 0; i < vec_size; i += inc)
        {
            b_type p1_vec = b_type::load_unaligned(&p1[i]);
            b_type t1_vec = b_type::load_unaligned(&t1[i]);
            b_type p2_vec = b_type::load_unaligned(&p2[i]);
            b_type t2_vec = b_type::load_unaligned(&t2[i]);

            b_type r_vec = m_r / cp * xsimd::log(p2_vec / p1_vec) / xsimd::log(t2_vec / t1_vec);
            r_vec.store_unaligned(&res[i]);
        }

        for (std::size_t i = vec_size; i < size; ++i)
            res[i] = eff_poly(p1[i], t1[i], p2[i], t2[i]);
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::h(const T& t) const
    {
        return m_cp * t;
    }

    template <class D>
    template <class T>
    void IdealGas<D>::h(T* res, T* temperature, std::size_t size) const
    {
        using b_type = xsimd::batch<double, xsimd::avx2>;
        std::size_t inc = b_type::size;

        std::size_t vec_size = size - size % inc;
        b_type cp = m_cp;

        for (std::size_t i = 0; i < vec_size; i += inc)
        {
            b_type t_vec = b_type::load_unaligned(&temperature[i]);
            b_type r_vec = cp * t_vec;
            r_vec.store_unaligned(&res[i]);
        }

        for (std::size_t i = vec_size; i < size; ++i)
            res[i] = m_cp * temperature[i];

        // thread_pool<std::size_t> pool(4);
        // // pool.resume();

        // auto vectorized_compute = [temperature, res, this](std::size_t /*runner*/, std::size_t start, std::size_t
        // end)
        // {
        //     using b_type = xsimd::batch<double, xsimd::avx2>;
        //     std::size_t inc = b_type::size;

        //     std::size_t size = end - start;
        //     std::size_t vec_end = start + size - size % inc;

        //     // std::cout << start << " / " << end << " / " << size << " / " << vec_end << std::endl << std::flush;
        //     for (std::size_t i = start; i < vec_end; i += inc)
        //     {
        //         b_type t_vec = b_type::load_unaligned(&temperature[i]);
        //         b_type r_vec = m_cp * t_vec;
        //         r_vec.store_unaligned(&res[i]);
        //     }

        //     for (std::size_t i = vec_end; i < end; ++i)
        //         res[i] = m_cp * temperature[i];
        // };

        // pool.run_blocks(0, size, vectorized_compute, 4);

        // pool.pause();
    }

    template <class D>
    inline double IdealGas<D>::r() const
    {
        return m_r;
    }

    template <class D>
    inline double IdealGas<D>::static_t(const double tt, const double mach) const
    {
        return tt / (1 + 0.5 * (m_cp / (m_cp - m_r) - 1.) * std::pow(mach, 2.));
    }

    template <class D>
    inline double IdealGas<D>::static_t(const double tt, const double mach, double, std::size_t) const
    {
        return static_t(tt, mach);
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::t_f_pr(const T& pr, const T& t1, const T& eff_poly) const
    {
        return t1 * std::pow(pr, m_r / (m_cp * eff_poly));
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::t_f_pr(const T& pr, const T& t1, const T& eff_poly, double, std::size_t) const
    {
        return t_f_pr(pr, t1, eff_poly);
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::t_f_h(const T& h) const
    {
        return h / m_cp;
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::t_f_h(const T& h, double, std::size_t) const
    {
        return t_f_h(h);
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::t_f_phi(const T& phi) const
    {
        return std::exp(phi / m_cp);
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::t_f_phi(const T& phi, double, std::size_t) const
    {
        return t_f_phi(phi);
    }

    template <class D>
    template <class T>
    auto IdealGas<D>::mach_f_wqa(const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter) const
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

    template <class D>
    inline bool operator==(const IdealGas<D>& ref, const IdealGas<D>& other)
    {
        return ref.m_cp == other.m_cp && ref.m_r == other.m_r;
    }
}

// clang-format off
#ifdef LIBTHERMO_USE_XTENSOR
    #include "libthermo/detail/ideal_gas_xt_impl.hpp"
#endif

#endif
