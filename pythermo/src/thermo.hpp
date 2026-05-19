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
    template <class A>
    class PyThermo
        : public thermo::ThermoExtendedInterface<double>
        , public thermo::ThermoInterface<A>
    {
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

    void thermo_base(nanobind::module_& m);

    void ideal_gas(nanobind::module_& m);

    void poly_gas(nanobind::module_& m);
}

#endif
