// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef PYTHERMO_GAS_H
#define PYTHERMO_GAS_H

#include "libthermo/thermo.hpp"
#include "libthermo/ideal_gas.hpp"
#include "libthermo/real_gas.hpp"

#include <pybind11/pybind11.h>
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


using array_t = xt::pyarray<double, xt::layout_type::row_major>;

class PyGas : public thermo::ThermoInterface<array_t>
{
public:
    using thermo::ThermoInterface<array_t>::ThermoInterface;

    double Gamma(double t) const override;
    array_t Gamma(const array_t& t) const override;

    double Cp(double t) const override;
    array_t Cp(const array_t& t);

    double H(double t) const override;
    array_t H(const array_t& t) const override;

    double Phi(double t) const override;
    array_t Phi(const array_t& t) const override;

    double PR(double t1, double t2, double eff_poly) const override;
    array_t PR(const array_t& t1, const array_t& t2, const array_t& eff_poly) const override;

    double EffPoly(double p1, double t1, double p2, double t2) const override;
    array_t EffPoly(const array_t& p1,
                    const array_t& t1,
                    const array_t& p2,
                    const array_t& t2) const override;

    double R() const override;

    double TFromPR(double pr_, double t1_, double eff_poly_) const override;

    double TFromH(double h_) const override;

    double TFromPhi(double phi_) const override;

    double StaticT(double tt_, double mach_) const override;
};

void
gas(pybind11::module_& m);

#endif
