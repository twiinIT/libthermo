// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"
#include "libthermo/ideal_gas.hpp"


namespace py = pybind11;
using namespace thermo;

using PyIdealGas = PyThermoHelper<IdealGas>;

void
ideal_gas(py::module_& m)
{
    using array_t = xt::pyarray<double, xt::layout_type::row_major>;
    using namespace py::literals;

    auto ideal_gas = py::class_<PyIdealGas, PyThermoInterface, std::shared_ptr<PyIdealGas>>(m, "IdealGas")
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

    /*
                  .def("static_t",
                       py::overload_cast<const double, const double, const double, const
       std::size_t>( &PyIdealGas::static_t, py::const_), "Static temperature", "tt"_a, "mach"_a,
                       "tol"_a,
                       "max_iter"_a = 30)
                  //.def("t_from_h", &PyIdealGas::t_f_h, "Temperature from enthalpy", "enthalpy"_a)
                  .def("t_from_h",
                       &PyIdealGas::t_f_h,
                       "Temperature from enthalpy",
                       "enthalpy"_a,
                       "tol"_a,
                       "max_iter"_a = 30)
                  //.def("t_from_phi", &PyIdealGas::t_f_phi, "Temperature from entropy",
       "entropy"_a) .def("t_from_phi", &PyIdealGas::t_f_phi, "Temperature from entropy",
                       "entropy"_a,
                       "tol"_a,
                       "max_iter"_a = 30)
                  //    .def("t_from_pr",
                  //         &PyIdealGas::t_f_pr,
                  //         "Temperature from pressure ratio, initial temperature and polytropic
                  //         efficiency", "pr"_a, "t1"_a, "eff_poly"_a)
                  .def("t_from_pr",
                       &PyIdealGas::t_f_pr,
                       "Temperature from pressure ratio, initial temperature and polytropic
       efficiency", "pr"_a, "t1"_a, "eff_poly"_a, "tol"_a, "max_iter"_a = 30);
    */
    // ideal_gas.def("static_t", &PyIdealGas::static_t, "Static temperature", "tt"_a, "mach"_a);
}
