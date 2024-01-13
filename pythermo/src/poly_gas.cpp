// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include "libthermo/poly_gas.hpp"


namespace py = pybind11;
using namespace thermo;

namespace pythermo
{
    using PG8 = PolyGas<PolyGasProps<double, 8>>;

    class PyPG8 : public PyThermoHelper<PG8, array_t>
    {
    public:
        using PG8::properties;

        PyPG8(const PolyGasProps<double, 8>& props)
            : PyThermoHelper<PG8, array_t>(props)
        {
        }
    };

    void poly_gas(py::module_& m)
    {
        using namespace py::literals;

        m.def("mix8", &mix<double, 8>, "cps"_a, "h0s"_a, "rs"_a, "weights"_a);

        auto polygasprops8 = py::class_<PolyGasProps<double, 8>>(m, "PolyGasProps8")
                                 .def(py::init<const std::array<double, 8>&, const double&, const double&>(),
                                      "cp_coeffs"_a,
                                      "h0"_a,
                                      "r"_a);
        polygasprops8.def(py::pickle(
            [](const PolyGasProps<double, 8>& p) -> py::tuple {  // __getstate__
                return py::make_tuple(p.cp_coeffs, p.h_coeffs[8], p.r);
            },
            [](py::tuple t) -> PolyGasProps<double, 8>
            {  // __setstate__
                PolyGasProps<double, 8> p(t[0].cast<std::array<double, 8>>(), t[1].cast<double>(), t[2].cast<double>());
                return p;
            }));

        auto polygas8 = py::class_<PyPG8, PyThermo<array_t>, std::shared_ptr<PyPG8>>(m, "PolyGas8");
        polygas8.def(py::init<const PolyGasProps<double, 8>&>(), "gas_properties"_a);
        polygas8.def(py::pickle(
            [](const PyPG8& g) -> py::tuple {  // __getstate__
                return py::make_tuple(&g.properties());
            },
            [](py::tuple t) -> PyPG8 {  // __setstate__
                PyPG8 g(t[0].cast<PolyGasProps<double, 8>>());
                return g;
            }));
    }
}