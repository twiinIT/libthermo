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
    class PyIdealGas : public PyThermo<IdealGas<double>>
    {
    public:
        PyIdealGas(double r, double cp)
        {
            p_gas = std::make_unique<IdealGas<double>>(r, cp);
        }

        double static_t_no_iter(const double& tt, const double& mach) const
        {
            return p_gas->static_t(tt, mach);
        }

        double t_f_h_no_iter(const double& h) const
        {
            return p_gas->t_f_h(h);
        }

        double t_f_phi_no_iter(const double& phi) const
        {
            return p_gas->t_f_phi(phi);
        }

        double t_f_pr_no_iter(const double& pr, const double& t1, const double& eff_poly) const
        {
            return p_gas->t_f_pr(pr, t1, eff_poly);
        }
    };

    void bind_ideal_gas(nb::module_& m)
    {
        using namespace nb::literals;

        auto gas = nb::class_<PyIdealGas>(m, "IdealGas").def(nb::init<double, double>());

        bind_thermo_extended_interface<PyIdealGas>(gas);

        gas.def("__getstate__", [](const PyIdealGas& g) -> nb::tuple { return nb::make_tuple(g.r(), g.cp(288.15)); })
            .def("__setstate__",
                 [](PyIdealGas& g, const nb::tuple& state) -> void
                 { new (&g) PyIdealGas(nb::cast<double>(state[0]), nb::cast<double>(state[1])); });

        gas.def("static_t", &PyIdealGas::static_t_no_iter, "Static temperature", "tt"_a, "mach"_a);
        gas.def("t_f_h", &PyIdealGas::t_f_h_no_iter, "Temperature from enthalpy", "h"_a);
        gas.def("t_f_phi", &PyIdealGas::t_f_phi_no_iter, "Temperature from phi function", "h"_a);
        gas.def("t_f_pr", &PyIdealGas::t_f_pr_no_iter, "Temperature from pressure ratio", "pr"_a, "t1"_a, "eff_poly"_a);
    }
}