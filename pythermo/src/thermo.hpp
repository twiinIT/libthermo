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

#include <initializer_list>

// #include <nanobind/stl.h>
// #include <nanobind/stl_bind.h>

namespace pythermo
{
    template <typename Class, typename Func, typename... Extra>
    NB_INLINE Class& multi_def(Class& class_, std::initializer_list<const char*> names, Func&& f, const Extra&... extra)
    {
        for (auto name : names)
        {
            class_.def(name, std::forward<Func>(f), extra...);
        }
        return class_;
    }

    template <typename Class, class ValueType>
    void register_base_interface(Class& class_)
    {
        using G = typename Class::Type;
        using namespace nanobind::literals;

        multi_def(class_, { "constant", "r" }, &G::r, "Gas constant");
        multi_def(class_, { "enthalpy", "h" }, &G::template h<const ValueType&>, "Enthalpy", "temperature"_a);
        multi_def(class_, { "entropy", "phi" }, &G::template phi<ValueType>, "Entropy", "temperature"_a);
        multi_def(class_,
                  { "specific_heat_ratio", "gamma" },
                  &G::template gamma<ValueType>,
                  "Specific heat ratio",
                  "temperature"_a);
        multi_def(class_,
                  { "specific_heat_pressure", "cp" },
                  &G::template cp<ValueType>,
                  "Specific heat pressure",
                  "temperature"_a);
        multi_def(class_,
                  { "pressure_ratio", "pr" },
                  &G::template pr<ValueType>,
                  "Pressure Ratio",
                  "t1"_a,
                  "t2"_a,
                  "eff_poly"_a);
        multi_def(class_,
                  { "polytropic_efficiency", "eff_poly" },
                  &G::template eff_poly<const ValueType&>,
                  "Polytropic efficiency",
                  "p1"_a,
                  "t1"_a,
                  "p2"_a,
                  "t2"_a);
    }

    template <typename Class, class ValueType>
    void register_extended_interface(Class& class_)
    {
        // def("static_t",
        //       &ThermoExtendedInterface<double>::static_t,
        //       "Static temperature",
        //       "tt"_a,
        //       "mach"_a,
        //       "tol"_a,
        //       "max_iter"_a = 30);

        // def("t_f_h",
        //       &ThermoExtendedInterface<double>::t_f_h,
        //       "Temperature from enthalpy",
        //       "h"_a,
        //       "tol"_a,
        //       "max_iter"_a = 30);

        // def("t_f_phi",
        //       &ThermoExtendedInterface<double>::t_f_phi,
        //       "Temperature from phi function",
        //       "h"_a,
        //       "tol"_a,
        //       "max_iter"_a = 30);

        // def("t_f_pr",
        //       &ThermoExtendedInterface<double>::t_f_pr,
        //       "Temperature from pressure ratio",
        //       "pr"_a,
        //       "t1"_a,
        //       "eff_poly"_a,
        //       "tol"_a,
        //       "max_iter"_a = 30);

        // def("mach_f_wqa",
        //       &ThermoExtendedInterface<double>::mach_f_wqa,
        //       "Mach number from specific mass flow",
        //       "pt"_a,
        //       "tt"_a,
        //       "wqa"_a,
        //       "tol"_a,
        //       "max_iter"_a = 30);
    }

    void thermo_base(nanobind::module_& m);

    void ideal_gas(nanobind::module_& m);

    void poly_gas(nanobind::module_& m);
}

#endif
