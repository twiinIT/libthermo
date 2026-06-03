// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include <nanobind/nanobind.h>

// using namespace thermo;
using namespace pythermo;

NB_MODULE(pythermo_core, m)
{
    m.doc() = "Python bindings of the C++ implementation of the thermodynamic library 'libthermo'.";

    bind_thermo_base(m);
    bind_ideal_gas(m);
    bind_poly_gas(m);
}
