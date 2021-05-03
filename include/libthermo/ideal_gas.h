// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_IDEAL_GAS_H
#define LIBTHERMO_IDEAL_GAS_H

#include "libthermo/gas.h"

#include <string>


namespace libthermo
{
    class IdealGas : public Gas<IdealGas>
    {
    public:
        IdealGas(double r_, double cp_);

        std::string Name() const;

        template<class T>
        T Gamma(const T& t) const;

        template<class T>
        T Cp(const T& = 0.) const;

        template<class T>
        T H(const T& t) const;
        
        template<class T>
        T Phi(const T&t) const;
        
        template<class T>
        T R() const;

        template<class T>
        T PR(const T& t1, const T& t2, const T& eff_poly) const;
        
        template<class T>
        T EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const;
        
        template<class T>
        T TFromH(const T& h) const;

    protected:
        double r, cp, gamma;
    };

    template<class T>
    T IdealGas::Gamma(const T& t) const
    {
        return gamma;
    }

    template<class T>
    T IdealGas::Cp(const T&) const
    {
        return cp;
    }
    
    template<class T>
    T IdealGas::H(const T& t) const
    {
        return cp * t;
    }

    template<class T>
    T IdealGas::Phi(const T& t) const
    {
        return cp * std::log(t);
    }

    template<class T>
    T IdealGas::R() const
    {
        return r; 
    }

    template<class T>
    T IdealGas::PR(const T& t1, const T& t2, const T& eff_poly) const
    { 
        return std::exp(std::log(t2 / t1) * eff_poly * cp / r);
    }

    template<class T>
    T IdealGas::EffPoly(const T& p1, const T& t1, const T& p2, const T& t2) const
    {
        return R<T>() * log(p2 / p1) / (Phi(t2) - Phi(t1));
    }
    
    template<class T>
    T IdealGas::TFromH(const T& h) const
    { 
        return h / cp; 
    }
}
#endif
