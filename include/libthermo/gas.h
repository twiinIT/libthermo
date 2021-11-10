// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_GAS_H
#define LIBTHERMO_GAS_H

#ifdef LIBTHERMO_USE_XTENSOR
#include "xtensor/xcontainer.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xfunction.hpp"

#define IS_NOT_XTENSOR \
    std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0
#define IS_XTENSOR \
    std::enable_if_t<xt::detail::is_container<T>::value, int> = 0
#else
#include <type_traits>

namespace
{
    template<class T>
    struct is_tensor : std::false_type {};

    template<class T>
    struct is_not_tensor : std::true_type {};
}
#define IS_NOT_XTENSOR \
    std::enable_if_t<is_not_tensor<T>::value, int> = 0
#define IS_XTENSOR \
    std::enable_if_t<is_tensor<T>::value, int> = 0
#endif

#include <cmath>
#include <string>
#include <memory>


namespace libthermo
{
    template <class G>
    class Gas
    {
    public:

        using derived_gas_type = G;

        inline const std::string& name() const { return m_name; };

        /// Proxy method to compute the specific heat ratio.
        template<class T>
        auto SpecificHeatRatio(const T& t) const { return derived_gas().Gamma(t); };

        /// Proxy method to compute the specific heat pressure.
        template<class T>
        auto SpecificHeatPressure(const T& t) const { return derived_gas().Cp(t); };

        /// Proxy method to compute the enthalpy.
        template<class T>
        auto Enthalpy(double const &t) const { return derived_gas().H(t); };
        
        /// Proxy method to compute the entropy.
        template<class T>
        auto Entropy(double const &t) const { return derived_gas().Phi(t); };
        
        /// Proxy method to get the gas constant.
        template<class T>
        auto Constant() const { return derived_gas().template R<T>(); };
        
        /// Proxy method to get the pressure ratio from initiale and finale temps, and polytropic efficiency.
        template<class T, class E = T>
        auto PressureRatio(const T& t1, const T& t2, const E& eff_poly) const { return derived_gas().PR(t1, t2, eff_poly); };
        
        /// Proxy method to get the polytropic efficiency of a transformation between initial and final pressure and temps.
        template<class T>
        auto PolytropicEfficiency(const T& p1, const T& t1, const T& p2, const T& t2) const { 
            return derived_gas().EffPoly(p1, t1, p2, t2); } 

    protected:

        Gas() = default;
        Gas(const std::string& name_);
        ~Gas() = default;

        const derived_gas_type& derived_gas() const noexcept;
        derived_gas_type& derived_gas() noexcept;

        const std::string m_name = "Gas";
    };

    template <class G>
    inline Gas<G>::Gas(const std::string& name_)
    : m_name(name_)
    {
    }

    template <class G>
    inline auto Gas<G>::derived_gas() const noexcept
        -> const derived_gas_type&
    {
        return *static_cast<const derived_gas_type*>(this);
    }

    template <class G>
    inline auto Gas<G>::derived_gas() noexcept
        -> derived_gas_type&
    {
        return *static_cast<derived_gas_type*>(this);
    }
}

#endif
