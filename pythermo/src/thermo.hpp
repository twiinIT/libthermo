// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef PYTHERMO_THERMO_H
#define PYTHERMO_THERMO_H

#include "libthermo/thermo.hpp"
#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"

#include <nanobind/nanobind.h>

// #include <nanobind/stl.h>
// #include <nanobind/stl_bind.h>

namespace pythermo
{
    template <class T>
    class PyThermo : public thermo::ThermoExtendedInterface<double>
    {
    protected:
        PyThermo() = default;
        std::unique_ptr<T> p_gas;

    public:
        virtual double r() const override
        {
            return p_gas->r();
        }
        virtual double gamma(const double& t) const override
        {
            return p_gas->gamma(t);
        }
        virtual double cp(const double& t) const override
        {
            return p_gas->cp(t);
        }
        virtual double phi(const double& t) const override
        {
            return p_gas->phi(t);
        }
        virtual double pr(const double& t1, const double& t2, const double& eff_poly) const override
        {
            return p_gas->pr(t1, t2, eff_poly);
        }
        virtual double eff_poly(const double& p1, const double& t1, const double& p2, const double& t2) const override
        {
            return p_gas->eff_poly(p1, t1, p2, t2);
        }
        virtual double h(const double& t) const override
        {
            return p_gas->h(t);
        }
        virtual double static_t(const double& tt,
                                const double& mach,
                                double tol,
                                std::size_t max_iter = 30) const override
        {
            return p_gas->static_t(tt, mach, tol, max_iter);
        }
        virtual double t_f_h(const double& h, double tol, std::size_t max_iter = 30) const override
        {
            return p_gas->t_f_h(h, tol, max_iter);
        }
        virtual double t_f_phi(const double& phi, double tol, std::size_t max_iter = 30) const override
        {
            return p_gas->t_f_phi(phi, tol, max_iter);
        }
        virtual double t_f_pr(const double& pr,
                              const double& t1,
                              const double& eff_poly,
                              double tol,
                              std::size_t max_iter = 30) const override
        {
            return p_gas->t_f_pr(pr, t1, eff_poly, tol, max_iter);
        }
        virtual double mach_f_wqa(
            const double& pt, const double& tt, const double& wqa, double tol, std::size_t max_iter = 30) const override
        {
            return p_gas->mach_f_wqa(pt, tt, wqa, tol, max_iter);
        }
    };

    template <class T, class C>
    void bind_thermo_extended_interface(nanobind::class_<C>& g)
    {
        using namespace nanobind::literals;
        using namespace thermo;

        g.def_prop_ro("constant", &T::r, "Gas constant");
        g.def("r", &T::r, "Gas constant");

        g.def("enthalpy", &T::h, "Enthalpy", "temperature"_a);
        g.def("h", &T::h, "Enthalpy", "temperature"_a);

        g.def("entropy", &T::phi, "Entropy", "temperature"_a);
        g.def("phi", &T::phi, "Entropy", "temperature"_a);

        g.def("specific_heat_ratio", &T::gamma, "Specific heat ratio", "temperature"_a);
        g.def("gamma", &T::gamma, "Specific heat ratio", "temperature"_a);

        g.def("specific_heat_pressure", &T::cp, "Specific heat pressure", "temperature"_a);
        g.def("cp", &T::cp, "Specific heat pressure", "temperature"_a);

        g.def("pressure_ratio", &T::pr, "Pressure Ratio", "t1"_a, "t2"_a, "eff_poly"_a);
        g.def("pr", &T::pr, "Pressure Ratio", "t1"_a, "t2"_a, "eff_poly"_a);

        g.def("polytropic_efficiency", &T::eff_poly, "Polytropic efficiency", "p1"_a, "t1"_a, "p2"_a, "t2"_a);
        g.def("eff_poly", &T::eff_poly, "Polytropic efficiency", "p1"_a, "t1"_a, "p2"_a, "t2"_a);

        g.def(
            "static_t",
            nanobind::overload_cast<const double&, const double&, double, std::size_t>(&T::static_t, nanobind::const_),
            "Static temperature",
            "tt"_a,
            "mach"_a,
            "tol"_a,
            "max_iter"_a = 30);

        g.def("t_f_h",
              nanobind::overload_cast<const double&, double, std::size_t>(&T::t_f_h, nanobind::const_),
              "Temperature from enthalpy",
              "h"_a,
              "tol"_a,
              "max_iter"_a = 30);

        g.def("t_f_phi",
              nanobind::overload_cast<const double&, double, std::size_t>(&T::t_f_phi, nanobind::const_),
              "Temperature from phi function",
              "h"_a,
              "tol"_a,
              "max_iter"_a = 30);

        g.def("t_f_pr",
              nanobind::overload_cast<const double&, const double&, const double&, double, std::size_t>(
                  &T::t_f_pr, nanobind::const_),
              "Temperature from pressure ratio",
              "pr"_a,
              "t1"_a,
              "eff_poly"_a,
              "tol"_a,
              "max_iter"_a = 30);

        g.def("mach_f_wqa",
              &T::mach_f_wqa,
              "Mach number from specific mass flow",
              "pt"_a,
              "tt"_a,
              "wqa"_a,
              "tol"_a,
              "max_iter"_a = 30);
    }

    void bind_thermo_base(nanobind::module_& m);

    void bind_ideal_gas(nanobind::module_& m);

    void bind_poly_gas(nanobind::module_& m);
}

#endif
