// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_THERMO_HPP
#define LIBTHERMO_THERMO_HPP

#include "libthermo/type_traits.hpp"

#include <cmath>
#include <string>
#include <memory>


namespace thermo
{
    template <class Th>
    class Thermo
    {
    public:
        using derived_type = Th;

        /// Proxy method to compute the specific heat ratio.
        template <class T>
        auto specific_heat_ratio(const T& t) const
        {
            return derived_thermo().gamma(t);
        };

        /// Proxy method to compute the specific heat pressure.
        template <class T>
        auto specific_heat_pressure(const T& t) const
        {
            return derived_thermo().cp(t);
        };

        /// Proxy method to compute the enthalpy.
        template <class T>
        auto enthalpy(const T& t) const
        {
            return derived_thermo().h(t);
        };

        /// Proxy method to compute the entropy.
        template <class T>
        auto entropy(const T& t) const
        {
            return derived_thermo().phi(t);
        };

        /// Proxy method to get the gas constant.
        double constant() const
        {
            return derived_thermo().r();
        };

        /// Proxy method to get the pressure ratio from initiale and finale temps, and polytropic
        /// efficiency.
        template <class T, class E = T>
        auto pressure_ratio(const T& t1, const T& t2, const E& eff_poly) const
        {
            return derived_thermo().pr(t1, t2, eff_poly);
        };

        /// Proxy method to get the polytropic efficiency of a transformation between initial and
        /// final pressure and temps.
        template <class T>
        auto polytropic_efficiency(const T& p1, const T& t1, const T& p2, const T& t2) const
        {
            return derived_thermo().eff_poly(p1, t1, p2, t2);
        }

    protected:
        Thermo() = default;
        ~Thermo() = default;

        const derived_type& derived_thermo() const noexcept;
        derived_type& derived_thermo() noexcept;
    };


    template <class T>
    inline auto Thermo<T>::derived_thermo() const noexcept -> const derived_type&
    {
        return *static_cast<const derived_type*>(this);
    }

    template <class G>
    inline auto Thermo<G>::derived_thermo() noexcept -> derived_type&
    {
        return *static_cast<derived_type*>(this);
    }


    template <class Vec>
    class ThermoInterface : public Thermo<ThermoInterface<Vec>>
    {
    public:
        ThermoInterface() = default;

        virtual double gamma(double t) const = 0;
        virtual Vec gamma(const Vec& t) const = 0;

        virtual double cp(double t) const = 0;
        virtual Vec cp(const Vec& t) const = 0;

        virtual double phi(double t) const = 0;
        virtual Vec phi(const Vec& t) const = 0;

        virtual double pr(double t1, double t2, double eff_poly) const = 0;
        virtual Vec pr(const Vec& t1, const Vec& t2, const Vec& eff_poly) const = 0;

        virtual double eff_poly(double p1, double t1, double p2, double t2) const = 0;
        virtual Vec eff_poly(const Vec& p1, const Vec& t1, const Vec& p2, const Vec& t2) const = 0;

        virtual double h(double t) const = 0;
        virtual Vec h(const Vec& t) const = 0;

        virtual double static_t(double tt, double mach, double tol, std::size_t max_iter) const = 0;
        // virtual Vec static_t(const Vec& tt, const Vec& mach) const = 0;

        virtual double t_f_h(double h, double tol, std::size_t max_iter) const = 0;
        // virtual Vec t_f_h(const Vec& h) const = 0;

        virtual double t_f_phi(double phi, double tol, std::size_t max_iter) const = 0;
        //  virtual Vec t_f_phi(const Vec& phi) const = 0;

        virtual double t_f_pr(
            double pr, double t1, double eff_poly, double tol, std::size_t max_iter) const = 0;
        // virtual Vec t_f_pr(const Vec& pr, const Vec& t1, const Vec& eff_poly) const = 0;

        virtual double mach_f_wqa(
            double pt, double tt, double wqa, double tol, std::size_t max_iter) const = 0;

        virtual double r() const = 0;
    };
}

#endif
