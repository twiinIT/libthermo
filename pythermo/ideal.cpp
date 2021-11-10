// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "libthermo/ideal_gas.h"
#include "libthermo/real_gas.h"

#include <pybind11/pybind11.h>
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <iostream>
#include <chrono>
#include <functional>


namespace py = pybind11;
using namespace libthermo;


PYBIND11_MODULE(pythermo, m)
{
    m.doc() = "Python bindings of the C++ implementation of the thermodynamic library 'libthermo'.";

    xt::import_numpy();

    using array_type = xt::pyarray<double, xt::layout_type::row_major>;

    py::class_<IdealGas>(m, "IdealGas")
        .def(py::init<double, double>())
        .def("constant", &IdealGas::Constant<double>, "Gas constant")
        .def("r", &IdealGas::R<double>, "Gas constant")
        .def("enthalpy", &IdealGas::Enthalpy<double>, "Enthalpy", py::arg("temperature"))
        .def("h", &IdealGas::H<double>, "Enthalpy", py::arg("temperature"))
        .def("entropy", &IdealGas::Entropy<double>, "Entropy", py::arg("temperature"))
        .def("phi", &IdealGas::Phi<double>, "Entropy", py::arg("temperature"))
        .def("phi", &IdealGas::Phi<array_type>, "Entropy", py::arg("temperature"))
        //.def("specific_heat_ratio", &IdealGas::SpecificHeatRatio<double>, "Specific heat ratio",
        // py::arg("temperature")) .def("gamma", &IdealGas::Gamma<double>, "Specific heat ratio",
        // py::arg("temperature"))
        .def("specific_heat_pressure",
             &IdealGas::SpecificHeatPressure<double>,
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp", &IdealGas::Cp<double>, "Specific heat pressure", py::arg("temperature"))
        .def("cp", &IdealGas::Cp<array_type>, "Specific heat pressure", py::arg("temperature"))
        .def("pressure_ratio",
             &IdealGas::PressureRatio<double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        //.def("pr", &IdealGas::PR<double>, "Pressure Ratio", py::arg("initiale temperature"),
        // py::arg("finale temperature"), py::arg("polytropic efficiency"))
        .def("pr",
             &IdealGas::PR<array_type, array_type>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             &IdealGas::PR<array_type, double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("eff_poly",
             &IdealGas::EffPoly<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("eff_poly",
             &IdealGas::PolytropicEfficiency<array_type>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("polytropic_efficiency",
             &IdealGas::PolytropicEfficiency<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"));

    py::class_<RealGas>(m, "RealGas")
        .def(py::init<double>())
        .def("constant", &RealGas::Constant<double>, "Gas constant")
        .def("r", &RealGas::R<double>, "Gas constant")
        .def("enthalpy", &RealGas::Enthalpy<double>, "Enthalpy", py::arg("temperature"))
        .def("h", &RealGas::H<double>, "Enthalpy", py::arg("temperature"))
        .def(
            "h",
            [](const RealGas& self, const array_type& t1) -> array_type { return self.H(t1); },
            "Enthalpy",
            py::arg("temperature"))
        .def("entropy", &RealGas::Entropy<double>, "Entropy", py::arg("temperature"))
        .def("phi", &RealGas::Phi<double>, "Entropy", py::arg("temperature"))
        .def(
            "phi",
            [](const RealGas& self, const array_type& t1) -> array_type { return self.Phi(t1); },
            "Entropy",
            py::arg("temperature"))
        .def(
            "dphi",
            [](const RealGas& self, const array_type& t1, const array_type& t2) -> array_type {
                return self.dPhi(t1, t2);
            },
            "Entropy delta",
            py::arg("initiale temperature"),
            py::arg("final temperature"))
        .def("specific_heat_ratio",
             &RealGas::SpecificHeatRatio<double>,
             "Specific heat ratio",
             py::arg("temperature"))
        .def(
            "gamma",
            [](const RealGas& self, const array_type& t1) -> array_type { return self.Gamma(t1); },
            "Specific heat ratio",
            py::arg("temperature"))
        .def("gamma", &RealGas::Gamma<double>, "Specific heat ratio", py::arg("temperature"))
        .def("specific_heat_pressure",
             &RealGas::SpecificHeatPressure<double>,
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp", &RealGas::Cp<double>, "Specific heat pressure", py::arg("temperature"))
        .def(
            "cp",
            [](const RealGas& self, const array_type& t) -> array_type { return self.Cp(t); },
            "Specific heat pressure",
            py::arg("temperature"))
        .def("pressure_ratio",
             &RealGas::PressureRatio<double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        //.def("pr", &RealGas::PR<double>, "Pressure Ratio", py::arg("initiale temperature"),
        // py::arg("finale temperature"), py::arg("polytropic efficiency"))
        .def(
            "pr",
            [](const RealGas& self,
               const array_type& t1,
               const array_type& t2,
               const array_type& eff_poly) -> array_type { return self.PR(t1, t2, eff_poly); },
            "Pressure Ratio",
            py::arg("initiale temperature"),
            py::arg("finale temperature"),
            py::arg("polytropic efficiency"))
        .def("pr",
             &RealGas::PR<array_type, double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("eff_poly",
             &RealGas::EffPoly<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def(
            "eff_poly",
            [](const RealGas& self,
               const array_type& p1,
               const array_type& t1,
               const array_type& p2,
               const array_type& t2) -> array_type { return self.EffPoly(p1, t1, p2, t2); },
            "Polytropic efficiency",
            py::arg("initial pressure"),
            py::arg("initial temperature"),
            py::arg("final pressure"),
            py::arg("final temperature"))
        .def("polytropic_efficiency",
             &RealGas::PolytropicEfficiency<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"));
}
