// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include "libthermo/poly_gas.hpp"

#include <nanobind/stl/vector.h>
#include <nanobind/stl/array.h>


namespace nb = nanobind;
using namespace thermo;

namespace pythermo
{
    class PyPoly8Gas : public PyThermo<PolyGas<PolyGasProps<double, 8>>>
    {
    public:
        PyPoly8Gas(const PolyGasProps<double, 8>& properties)
        {
            p_gas = std::make_unique<PolyGas<PolyGasProps<double, 8>>>(properties);
        }

        const PolyGasProps<double, 8>& properties() const
        {
            return p_gas->properties();
        }
    };


    void bind_poly_gas(nb::module_& m)
    {
        using namespace nb::literals;

        auto polygasprops8 = nb::class_<PolyGasProps<double, 8>>(m, "PolyGasProps8")
                                 .def(nb::init<const std::array<double, 8>&, const double&, const double&>(),
                                      "cp_coeffs"_a,
                                      "h0"_a,
                                      "r"_a);

        polygasprops8
            .def("__getstate__",
                 [](const PolyGasProps<double, 8>& p) -> nb::tuple
                 { return nb::make_tuple(p.cp_coeffs, p.h_coeffs[8], p.r); })
            .def("__setstate__",
                 [](PolyGasProps<double, 8>& p, const nb::tuple& state) -> void
                 {
                     new (&p) PolyGasProps<double, 8>(nb::cast<std::array<double, 8>>(state[0]),
                                                      nb::cast<double>(state[1]),
                                                      nb::cast<double>(state[2]));
                 });

        m.def("mix8", &mix<double, 8>, "cps"_a, "h0s"_a, "rs"_a, "weights"_a);

        auto gas = nb::class_<PyPoly8Gas>(m, "PolyGas8");
        gas.def(nb::init<const PolyGasProps<double, 8>&>(), "gas_properties"_a);

        bind_thermo_extended_interface<PyPoly8Gas>(gas);

        gas.def("__getstate__", [](const PyPoly8Gas& g) -> nb::tuple { return nb::make_tuple(g.properties()); })
            .def("__setstate__",
                 [](PyPoly8Gas& g, const nb::tuple& state) -> void
                 { new (&g) PyPoly8Gas(nb::cast<PolyGasProps<double, 8>>(state[0])); });
    }
}