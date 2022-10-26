// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#define FORCE_IMPORT_ARRAY

#include "thermo.hpp"

#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

namespace py = pybind11;
using namespace thermo;
using array_type = xt::pyarray<double, xt::layout_type::row_major>;

/*
void
ideal_gas(py::module_&);
void
poly_thermo(py::module_&);
*/

/*
void
print_h(ThermoInterface<array_type>* gas)
{
    double h = gas->h(288.15);
    std::cout << "h: " << h << std::endl;
}
*/

PYBIND11_MODULE(pythermo, m)
{
    xt::import_numpy();

    m.doc() = "Python bindings of the C++ implementation of the thermodynamic library 'libthermo'.";
    thermo_base(m);
    ideal_gas(m);
    poly_gas(m);
}
