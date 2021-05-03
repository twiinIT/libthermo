// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMOCORE_GAS_H
#define LIBTHERMOCORE_GAS_H

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

        inline const std::string Name() const { return "Gas"; };

        /// Proxy method to compute the specific heat ratio.
        double SpecificHeatRatio(double const &t) const { return derived_gas().Gamma(t); };

        /// Proxy method to compute the specific heat pressure.
        double SpecificHeatPressure(double const &t) const { return derived_gas().Cp(t); };

        /// Proxy method to compute the enthalpy.
        double Enthalpy(double const &t) const { return derived_gas().H(t); };
        
        /// Proxy method to compute the entropy.
        double Entropy(double const &t) const { return derived_gas().Phi(t); };
        
        /// Proxy method to get the gas constant.
        double Constant() const { return derived_gas().R(); };
        
        /// Proxy method to get the pressure ratio from initiale and finale temps, and polytropic efficiency.
        double PressureRatio(const double &t1, const double &t2, const double &eff_poly) const { return derived_gas().PR(t1, t2, eff_poly); };
        
        /// Proxy method to get the polytropic efficiency of a transformation between initial and final pressure and temps.
        double PolytropicEfficiency(const double &p1, const double &t1, const double &p2, const double &t2) const { 
            return derived_gas().EffPoly(p1, t1, p2, t2); } 

    protected:

        Gas() = default;
        ~Gas() = default;

        const derived_gas_type& derived_gas() const noexcept;
        derived_gas_type& derived_gas() noexcept;
    };

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
