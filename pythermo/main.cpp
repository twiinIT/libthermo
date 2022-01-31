// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "gas.hpp"
#include "libthermo/ideal_gas.hpp"
#include "libthermo/real_gas.hpp"


namespace py = pybind11;
using namespace thermo;
using array_type = xt::pyarray<double, xt::layout_type::row_major>;

/*
void
ideal_gas(py::module_&);
void
real_gas(py::module_&);
*/

/*
void
print_h(ThermoInterface<array_type>* gas)
{
    double h = gas->H(288.15);
    std::cout << "h: " << h << std::endl;
}
*/

PYBIND11_MODULE(pythermo, m)
{
    m.doc() = "Python bindings of the C++ implementation of the thermodynamic library 'libthermo'.";

    xt::import_numpy();
    gas(m);
    // ideal_gas(m);
    // real_gas(m);

    // m.def("print_h", &print_h);
}
