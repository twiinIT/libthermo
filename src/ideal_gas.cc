// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "libthermo/ideal_gas.h"


namespace libthermo
{
    IdealGas::IdealGas(double r_, double cp_)
        : r(r_), cp(cp_), gamma(cp_ / (cp_ - r_)) {}

    std::string IdealGas::Name() const
    {
        return "IdealGas";
    }
}
