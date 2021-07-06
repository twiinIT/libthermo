// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "libthermo/ideal_gas.h"

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

  py::class_<IdealGas> (m, "IdealGas")
  .def(py::init<double, double>())
  .def("constant", &IdealGas::Constant<double>, "Gas constant")
  .def("r", &IdealGas::R<double>, "Gas constant")
  .def("enthalpy", &IdealGas::Enthalpy<double>, "Enthalpy", py::arg("temperature"))
  .def("h", &IdealGas::H<double>, "Enthalpy", py::arg("temperature"))
  .def("entropy", &IdealGas::Entropy<double>, "Entropy", py::arg("temperature"))
  .def("phi", &IdealGas::Phi<double>, "Entropy", py::arg("temperature"))
  .def("phi", &IdealGas::Phi<array_type>, "Entropy", py::arg("temperature"))
  //.def("specific_heat_ratio", &IdealGas::SpecificHeatRatio<double>, "Specific heat ratio", py::arg("temperature"))
  //.def("gamma", &IdealGas::Gamma<double>, "Specific heat ratio", py::arg("temperature"))
  .def("specific_heat_pressure", &IdealGas::SpecificHeatPressure<double>, "Specific heat pressure", py::arg("temperature"))
  .def("cp", &IdealGas::Cp<double>, "Specific heat pressure", py::arg("temperature"))
  .def("pressure_ratio", &IdealGas::PressureRatio<double>, "Pressure Ratio", py::arg("initiale temperature"), py::arg("finale temperature"), py::arg("polytropic efficiency"))
  //.def("pr", &IdealGas::PR<double>, "Pressure Ratio", py::arg("initiale temperature"), py::arg("finale temperature"), py::arg("polytropic efficiency"))
  .def("pr", &IdealGas::PR<array_type, array_type>, "Pressure Ratio", py::arg("initiale temperature"), py::arg("finale temperature"), py::arg("polytropic efficiency"))
  .def("pr", &IdealGas::PR<array_type, double>, "Pressure Ratio", py::arg("initiale temperature"), py::arg("finale temperature"), py::arg("polytropic efficiency"))
  .def("eff_poly", &IdealGas::EffPoly<double>, "Polytropic efficiency", 
        py::arg("initial pressure"), py::arg("initial temperature"), py::arg("final pressure"), py::arg("final temperature"))
  .def("eff_poly", &IdealGas::PolytropicEfficiency<array_type>, "Polytropic efficiency", 
        py::arg("initial pressure"), py::arg("initial temperature"), py::arg("final pressure"), py::arg("final temperature"))
  .def("polytropic_efficiency", &IdealGas::PolytropicEfficiency<double>, "Polytropic efficiency", 
        py::arg("initial pressure"), py::arg("initial temperature"), py::arg("final pressure"), py::arg("final temperature"));

}

