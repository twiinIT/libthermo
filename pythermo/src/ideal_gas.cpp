// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"
#include "libthermo/ideal_gas.hpp"


namespace py = pybind11;
using namespace thermo;

namespace pythermo
{
    class PyIdealGas : public PyThermoHelper<IdealGas<double>, array_t>
    {
    public:
        PyIdealGas(double r, double cp)
            : PyThermoHelper<IdealGas<double>, array_t>(r, cp)
        {
        }

        double static_t(const double& t, const double& mach) const
        {
            return IdealGas<double>::static_t(t, mach);
        }
        double t_f_h(const double& h) const
        {
            return IdealGas<double>::t_f_h(h);
        }
        double t_f_phi(const double& phi) const
        {
            return IdealGas<double>::t_f_phi(phi);
        }
        double t_f_pr(const double& pr, const double& t1, const double& eff_poly) const
        {
            return IdealGas<double>::t_f_pr(pr, t1, eff_poly);
        }
    };

    void ideal_gas(py::module_& m)
    {
        using namespace py::literals;

        auto ideal_gas = py::class_<PyIdealGas, PyThermo<array_t>, std::shared_ptr<PyIdealGas>>(m, "IdealGas")
                             .def(py::init<double, double>())
                             .def_property_readonly("constant", &PyIdealGas::r, "Gas constant")
                             .def_property_readonly("r", &PyIdealGas::r, "Gas constant");

        ideal_gas.def(py::pickle(
            [](const PyIdealGas& g) -> py::tuple {  // __getstate__
                return py::make_tuple(g.r(), g.cp(288.15));
            },
            [](py::tuple t) -> PyIdealGas {  // __setstate__
                PyIdealGas g(t[0].cast<double>(), t[1].cast<double>());
                return g;
            }));

        // Define `IdealGas` overloads of inverse methods not requiring iterations
        // Those overloads are NOT part of the `ThermoExtendedInterface` but are provided for convenience
        // note: base class overloads need to be redefined here to not be masked by derived class ones
        ideal_gas.def(
            "static_t", &PyIdealGas::static_t, "Static temperature from total and Mach number", "tt"_a, "mach"_a);
        ideal_gas.def("static_t",
                      &ThermoExtendedInterface<double>::static_t,
                      "Static temperature from total and Mach number",
                      "tt"_a,
                      "mach"_a,
                      "tol"_a,
                      "max_iter"_a = 30);

        ideal_gas.def("t_f_h", &PyIdealGas::t_f_h, "Temperature from enthalpy", "h"_a);
        ideal_gas.def("t_f_h",
                      &ThermoExtendedInterface<double>::t_f_h,
                      "Temperature from enthalpy",
                      "h"_a,
                      "tol"_a,
                      "max_iter"_a = 30);

        ideal_gas.def("t_f_phi", &PyIdealGas::t_f_phi, "Temperature from phi function", "phi"_a);
        ideal_gas.def("t_f_phi",
                      &ThermoExtendedInterface<double>::t_f_phi,
                      "Temperature from phi function",
                      "h"_a,
                      "tol"_a,
                      "max_iter"_a = 30);

        ideal_gas.def("t_f_pr", &PyIdealGas::t_f_pr, "Temperature from pressure ratio", "pr"_a, "t1"_a, "eff_poly"_a);
        ideal_gas.def("t_f_pr",
                      &ThermoExtendedInterface<double>::t_f_pr,
                      "Temperature from pressure ratio",
                      "pr"_a,
                      "t1"_a,
                      "eff_poly"_a,
                      "tol"_a,
                      "max_iter"_a = 30);

        ideal_gas.def("t_f_pr", &PyIdealGas::t_f_pr, "Temperature from pressure ratio", "pr"_a, "t1"_a, "eff_poly"_a);
        ideal_gas.def("t_f_pr",
                      &ThermoExtendedInterface<double>::t_f_pr,
                      "Temperature from pressure ratio",
                      "pr"_a,
                      "t1"_a,
                      "eff_poly"_a,
                      "tol"_a,
                      "max_iter"_a = 30);
    }
}