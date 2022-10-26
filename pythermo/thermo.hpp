// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef PYTHERMO_GAS_H
#define PYTHERMO_GAS_H

#include "libthermo/thermo.hpp"
#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"

#include <pybind11/pybind11.h>

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>


using array_t = xt::pyarray<double, xt::layout_type::row_major>;

class PyThermo : public thermo::ThermoInterface<array_t>
{
public:
    using thermo::ThermoInterface<array_t>::ThermoInterface;

    double gamma(double t) const override;
    array_t gamma(const array_t& t) const override;

    double cp(double t) const override;
    array_t cp(const array_t& t) const;

    double h(double t) const override;
    array_t h(const array_t& t) const override;

    double phi(double t) const override;
    array_t phi(const array_t& t) const override;

    double pr(double t1, double t2, double eff_poly) const override;
    array_t pr(const array_t& t1, const array_t& t2, const array_t& eff_poly) const override;

    double eff_poly(double p1, double t1, double p2, double t2) const override;
    array_t eff_poly(const array_t& p1,
                     const array_t& t1,
                     const array_t& p2,
                     const array_t& t2) const override;

    double r() const override;

    double t_f_pr(double pr_,
                  double t1_,
                  double eff_poly_,
                  double tol,
                  std::size_t max_iter = 30) const override;

    double t_f_h(double h_, double tol, std::size_t max_iter = 30) const override;

    double t_f_phi(double phi_, double tol, std::size_t max_iter = 30) const override;

    double static_t(double tt_, double mach_, double tol, std::size_t max_iter = 30) const override;

    double mach_f_wqa(
        double pt, double tt, double wqa, double tol, std::size_t max_iter = 30) const override;
};

void
thermo_base(pybind11::module_& m);

void
ideal_gas(pybind11::module_& m);

void
poly_gas(pybind11::module_& m);

#endif
