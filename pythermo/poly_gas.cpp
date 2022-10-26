// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include "libthermo/poly_gas.hpp"


namespace py = pybind11;
using namespace thermo;

void
poly_gas(py::module_& m)
{
    using array_t = xt::pyarray<double, xt::layout_type::row_major>;
    using PG8 = PolyGas<PolyGasProps<8>>;
    using namespace py::literals;

    py::class_<PolyGasProps<8>>(m, "PolyGasProps8")
        .def(
            py::init<const std::array<double, 8>&, double, double>(), "cp_coeffs"_a, "h0"_a, "r"_a);
    m.def("mix8", &mix<8>, "cps"_a, "h0s"_a, "rs"_a, "weights"_a);

    auto polygas8 = py::class_<PG8>(m, "PolyGas8");

    polygas8.def(py::init<const PolyGasProps<8>&>(), "gas_properties"_a)
        .def_property_readonly("constant", &PG8::r, "Gas constant")
        .def_property_readonly("r", &PG8::r, "Gas constant")
        .def("h", &PG8::h<double>, "Enthalpy", "temperature"_a)
        .def("gamma", &PG8::gamma<double>, "Specific heat ratio", "temperature"_a)
        .def("cp", &PG8::cp<double>, "Specific heat pressure", "temperature"_a)
        .def("phi", &PG8::phi<double>, "Entropy", "temperature"_a)
        .def("pr", &PG8::pr<double>, "Pressure Ratio", "t_in"_a, "t_out"_a, "eff_poly"_a)
        .def("eff_poly",
             &PG8::eff_poly<double>,
             "Polytropic efficiency",
             "p_in"_a,
             "p_out"_a,
             "t_in"_a,
             "t_out"_a)
        .def("static_t",
             &PG8::static_t<double>,
             "Static temperature",
             "tt"_a,
             "mach"_a,
             "tol"_a,
             "max_iter"_a = 30)
        .def("t_from_h",
             &PG8::t_f_h<double>,
             "Temperature from enthalpy",
             "h"_a,
             "tol"_a,
             "max_iter"_a = 30)
        .def("t_from_phi",
             &PG8::t_f_phi<double>,
             "Temperature from entropy",
             "s"_a,
             "tol"_a,
             "max_iter"_a = 30)
        .def("t_from_pr",
             &PG8::t_f_pr<double>,
             "Temperature from pressure ratio, initial temperature and polytropic efficiency",
             "pr"_a,
             "t1"_a,
             "eff_poly"_a,
             "tol"_a,
             "max_iter"_a = 30)
        .def("mach_f_wqa",
             &PG8::mach_f_wqa<double>,
             "Mach from mass flow over section ratio",
             "pt"_a,
             "tt"_a,
             "wqa"_a,
             "tol"_a,
             "max_iter"_a = 30);

    polygas8
        .def(
            "h",
            [](const PG8& gas, const array_t& t) -> array_t { return gas.h(t); },
            "Enthalpy",
            py::arg("temperature"))
        .def(
            "gamma",
            [](const PG8& gas, const array_t& t) -> array_t { return gas.gamma(t); },
            "Specific heat ratio",
            py::arg("temperature"))
        .def(
            "cp",
            [](const PG8& gas, const array_t& t) -> array_t { return gas.cp(t); },
            "Specific heat pressure",
            py::arg("temperature"))
        .def(
            "phi",
            [](const PG8& gas, const array_t& t) -> array_t { return gas.phi(t); },
            "Entropy",
            py::arg("temperature"))
        .def(
            "eff_poly",
            [](const PG8& gas,
               const array_t& pin,
               const array_t& pout,
               const array_t& tin,
               const array_t& tout) -> array_t { return gas.eff_poly(pin, pout, tin, tout); },
            "Polytropic efficiency",
            "p_in"_a,
            "p_out"_a,
            "t_in"_a,
            "t_out"_a)
        .def(
            "pr",
            [](const PG8& gas, const array_t& tin, const array_t& tout, const array_t& eff_poly)
                -> array_t { return gas.pr(tin, tout, eff_poly); },
            "Pressure Ratio",
            "t_in"_a,
            "t_out"_a,
            "eff_poly"_a);
}
