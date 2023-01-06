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

class PyThermoInterface
    : public thermo::ThermoExtendedInterface<double>
    , public thermo::ThermoInterface<array_t>
{
};


#define FORWARD_CALL(tmpl, name, type, ...)                                                                            \
    type name(const type& t) const override                                                                            \
    {                                                                                                                  \
        return tmpl::name(t);                                                                                          \
    }

template <class G>
class PyThermoHelper
    : public PyThermoInterface
    , protected G
{
public:
    using G::G;

    FORWARD_CALL(G, gamma, double);
    FORWARD_CALL(G, gamma, array_t);

    FORWARD_CALL(G, cp, double);
    FORWARD_CALL(G, cp, array_t);

    FORWARD_CALL(G, h, double);
    FORWARD_CALL(G, h, array_t);

    FORWARD_CALL(G, phi, double);
    FORWARD_CALL(G, phi, array_t);

    double r() const override
    {
        return G::r();
    }

    double pr(const double& t1, const double& t2, const double& eff_poly) const override
    {
        return G::pr(t1, t2, eff_poly);
    }
    array_t pr(const array_t& t1, const array_t& t2, const array_t& eff_poly) const override
    {
        return G::pr(t1, t2, eff_poly);
    }

    double eff_poly(const double& p1, const double& t1, const double& p2, const double& t2) const override
    {
        return G::eff_poly(p1, t1, p2, t2);
    }
    array_t eff_poly(const array_t& p1, const array_t& t1, const array_t& p2, const array_t& t2) const override
    {
        return G::eff_poly(p1, t1, p2, t2);
    }

    double static_t(const double& t, const double& mach) const
    {
        return G::static_t(t, mach);
    }

    double static_t(const double& t, const double& mach, double tol, const std::size_t max_iter = 30) const override
    {
        return G::static_t(t, mach, tol, max_iter);
    }

    double t_f_pr(const double& pr,
                  const double& t1,
                  const double& eff_poly,
                  double tol,
                  std::size_t max_iter = 30) const override
    {
        return G::t_f_pr(pr, t1, eff_poly, tol, max_iter);
    }
    double t_f_h(const double& h, double tol, std::size_t max_iter = 30) const override
    {
        return G::t_f_h(h, tol, max_iter);
    }
    double t_f_phi(const double& phi, double tol, std::size_t max_iter = 30) const override
    {
        return G::t_f_phi(phi, tol, max_iter);
    }
    double mach_f_wqa(
        const double& pt, const double& tt, const double& wqa, double tol, std::size_t max_iter = 30) const override
    {
        return G::mach_f_wqa(pt, tt, wqa, tol, max_iter);
    }
};


void
thermo_base(pybind11::module_& m);

void
ideal_gas(pybind11::module_& m);

void
poly_gas(pybind11::module_& m);

#endif
