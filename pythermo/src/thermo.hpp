// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef PYTHERMO_THERMO_H
#define PYTHERMO_THERMO_H

#include "libthermo/thermo.hpp"
#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"

#include <pybind11/pybind11.h>

#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace pythermo
{
    using array_t = xt::pyarray<double, xt::layout_type::row_major>;

    template <class A>
    class PyThermo
        : public thermo::ThermoExtendedInterface<double>
        , public thermo::ThermoInterface<A>
    {
    };

#define FORWARD_TO_BASE(tmpl, name, type, ...)                                                                         \
    type name(const type& t) const override                                                                            \
    {                                                                                                                  \
        return tmpl::name(t);                                                                                          \
    }

    template <class G, class A>
    class PyThermoHelper
        : public PyThermo<A>
        , protected G
    {
    public:
        using G::G;

        FORWARD_TO_BASE(G, gamma, double);
        FORWARD_TO_BASE(G, gamma, A);

        FORWARD_TO_BASE(G, cp, double);
        FORWARD_TO_BASE(G, cp, A);

        FORWARD_TO_BASE(G, h, double);
        FORWARD_TO_BASE(G, h, A);

        FORWARD_TO_BASE(G, phi, double);
        FORWARD_TO_BASE(G, phi, A);

        double r() const override
        {
            return G::r();
        }

        double pr(const double& t1, const double& t2, const double& eff_poly) const override
        {
            return G::pr(t1, t2, eff_poly);
        }
        A pr(const A& t1, const A& t2, const A& eff_poly) const override
        {
            return G::pr(t1, t2, eff_poly);
        }

        double eff_poly(const double& p1, const double& t1, const double& p2, const double& t2) const override
        {
            return G::eff_poly(p1, t1, p2, t2);
        }
        A eff_poly(const A& p1, const A& t1, const A& p2, const A& t2) const override
        {
            return G::eff_poly(p1, t1, p2, t2);
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

#undef FORWARD_TO_BASE

    void thermo_base(pybind11::module_& m);

    void ideal_gas(pybind11::module_& m);

    void poly_gas(pybind11::module_& m);
}

#endif
