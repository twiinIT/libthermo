// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_THERMO_HPP
#define LIBTHERMO_THERMO_HPP

#ifdef LIBTHERMO_USE_XTENSOR
#include "xtensor/xcontainer.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xfunction.hpp"

#define IS_NOT_XTENSOR std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0
#define IS_XTENSOR std::enable_if_t<xt::detail::is_container<T>::value, int> = 0
#else
#include <type_traits>

namespace
{
    template <class T>
    struct is_tensor : std::false_type
    {
    };

    template <class T>
    struct is_not_tensor : std::true_type
    {
    };
}
#define IS_NOT_XTENSOR std::enable_if_t<is_not_tensor<T>::value, int> = 0
#define IS_XTENSOR std::enable_if_t<is_tensor<T>::value, int> = 0
#endif

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

        inline const std::string& name() const
        {
            return m_name;
        };

        /// Proxy method to compute the specific heat ratio.
        template <class T>
        auto SpecificHeatRatio(const T& t) const
        {
            return derived_thermo().Gamma(t);
        };

        /// Proxy method to compute the specific heat pressure.
        template <class T>
        auto SpecificHeatPressure(const T& t) const
        {
            return derived_thermo().Cp(t);
        };

        /// Proxy method to compute the enthalpy.
        template <class T>
        auto Enthalpy(const T& t) const
        {
            return derived_thermo().H(t);
        };

        /// Proxy method to compute the entropy.
        template <class T>
        auto Entropy(const T& t) const
        {
            return derived_thermo().Phi(t);
        };

        /// Proxy method to get the gas constant.
        double Constant() const
        {
            return derived_thermo().R();
        };

        /// Proxy method to get the pressure ratio from initiale and finale temps, and polytropic
        /// efficiency.
        template <class T, class E = T>
        auto PressureRatio(const T& t1, const T& t2, const E& eff_poly) const
        {
            return derived_thermo().PR(t1, t2, eff_poly);
        };

        /*
                /// Proxy method to get the temperature from pressure ratio, initial temperature,
           and
                /// polytropic efficiency.
                template <class T, class E = T>
                auto TempFromPR(const T& pr, const T& p2, const E& eff_poly) const
                {
                    return derived_thermo().Tau(p1, p2, eff_poly);
                };
        */
        /// Proxy method to get the polytropic efficiency of a transformation between initial and
        /// final pressure and temps.
        template <class T>
        auto PolytropicEfficiency(const T& p1, const T& t1, const T& p2, const T& t2) const
        {
            return derived_thermo().EffPoly(p1, t1, p2, t2);
        }

    protected:
        Thermo() = default;
        Thermo(const std::string& name_);
        ~Thermo() = default;

        const derived_type& derived_thermo() const noexcept;
        derived_type& derived_thermo() noexcept;

        const std::string m_name = "";
    };

    template <class T>
    inline Thermo<T>::Thermo(const std::string& name_)
        : m_name(name_)
    {
    }

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

        virtual double Gamma(double t) const = 0;
        virtual Vec Gamma(const Vec& t) const = 0;

        virtual double Cp(double t) const = 0;
        virtual Vec Cp(const Vec& t) const = 0;

        virtual double Phi(double t) const = 0;
        virtual Vec Phi(const Vec& t) const = 0;

        virtual double PR(double t1, double t2, double eff_poly) const = 0;
        virtual Vec PR(const Vec& t1, const Vec& t2, const Vec& eff_poly) const = 0;

        virtual double EffPoly(double p1, double t1, double p2, double t2) const = 0;
        virtual Vec EffPoly(const Vec& p1, const Vec& t1, const Vec& p2, const Vec& t2) const = 0;

        virtual double H(double t) const = 0;
        virtual Vec H(const Vec& t) const = 0;

        virtual double StaticT(double tt, double mach) const = 0;
        // virtual Vec StaticT(const Vec& tt, const Vec& mach) const = 0;

        virtual double TFromH(double h) const = 0;
        // virtual Vec TFromH(const Vec& h) const = 0;

        virtual double TFromPhi(double phi) const = 0;
        //  virtual Vec TFromPhi(const Vec& phi) const = 0;

        virtual double TFromPR(double pr, double t1, double eff_poly) const = 0;
        // virtual Vec TFromPR(const Vec& pr, const Vec& t1, const Vec& eff_poly) const = 0;

        virtual double R() const = 0;
    };
}

#endif
