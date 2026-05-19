// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"
#include "libthermo/ideal_gas.hpp"


namespace nb = nanobind;
using namespace thermo;

namespace pythermo
{
    void ideal_gas(nb::module_& m)
    {
        using namespace nb::literals;

        auto gas = nb::class_<IdealGas<double>>(m, "IdealGas");
        gas.def(nb::init<double, double>());

        bind_thermo_extended_interface<IdealGas<double>>(gas);

        // ideal_gas.def(nb::pickle(
        //     [](const IdealGas<double>& g) -> nb::tuple {  // __getstate__
        //         return nb::make_tuple(g.r(), g.cp(288.15));
        //     },
        //     [](nb::tuple t) -> IdealGas<double> {  // __setstate__
        //         IdealGas<double> g(t[0].cast<double>(), t[1].cast<double>());
        //         return g;
        //     }));

        gas.def("static_t",
                nanobind::overload_cast<const double&, const double&>(&IdealGas<double>::static_t, nanobind::const_),
                "Static temperature",
                "tt"_a,
                "mach"_a);

        gas.def("t_f_h",
                nanobind::overload_cast<const double&>(&IdealGas<double>::t_f_h, nanobind::const_),
                "Temperature from enthalpy",
                "h"_a);

        gas.def("t_f_phi",
                nanobind::overload_cast<const double&>(&IdealGas<double>::t_f_phi, nanobind::const_),
                "Temperature from phi function",
                "h"_a);

        gas.def("t_f_pr",
                nanobind::overload_cast<const double&, const double&, const double&>(&IdealGas<double>::t_f_pr,
                                                                                     nanobind::const_),
                "Temperature from pressure ratio",
                "pr"_a,
                "t1"_a,
                "eff_poly"_a);
    }
}