// Copyright (c) 2021-2026, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include "libthermo/exceptions.hpp"

#include <nanobind/trampoline.h>


namespace nb = nanobind;
using namespace thermo;

namespace pythermo
{
#define THERMO_INTERFACE_BASE_TRAMPOLINE()                                                                             \
    double r() const override                                                                                          \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(r);                                                                                           \
    }

#define THERMO_INTERFACE_TRAMPOLINE(T)                                                                                 \
    T gamma(const T& t) const override                                                                                 \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(gamma, t);                                                                                    \
    }                                                                                                                  \
    T cp(const T& t) const override                                                                                    \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(cp, t);                                                                                       \
    }                                                                                                                  \
    T h(const T& t) const override                                                                                     \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(h, t);                                                                                        \
    }                                                                                                                  \
    T phi(const T& t) const override                                                                                   \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(phi, t);                                                                                      \
    }                                                                                                                  \
    T pr(const T& t1, const T& t2, const T& eff_poly) const override                                                   \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(pr, t1, t2, eff_poly);                                                                        \
    }                                                                                                                  \
    T eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const override                                      \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(eff_poly, p1, t1, p2, t2);                                                                    \
    }

#define THERMO_EXT_INTERFACE_TRAMPOLINE(T)                                                                             \
    T t_f_h(const T& h_, double tol, std::size_t max_iter = 30) const override                                         \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(t_f_h, h_, tol, max_iter);                                                                    \
    }                                                                                                                  \
    T t_f_phi(const T& phi_, double tol, std::size_t max_iter = 30) const override                                     \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(t_f_phi, phi_, tol, max_iter);                                                                \
    }                                                                                                                  \
    T t_f_pr(const T& pr_, const T& t1_, const T& eff_poly_, double tol, std::size_t max_iter = 30) const override     \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(t_f_pr, pr_, t1_, eff_poly_, tol, max_iter);                                                  \
    }                                                                                                                  \
    T static_t(const T& tt_, const T& mach_, double tol, std::size_t max_iter = 30) const override                     \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(static_t, tt_, mach_, tol, max_iter);                                                         \
    }                                                                                                                  \
    T mach_f_wqa(const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter = 30) const override         \
    {                                                                                                                  \
        NB_OVERRIDE_PURE(mach_f_wqa, pt, tt, wqa, tol, max_iter);                                                      \
    }                                                                                                                  \
    THERMO_INTERFACE_TRAMPOLINE(T)

    /* Defines the C++ trampolines to allow Python users to
     * define their own thermo and use it in code accepting
     * a PyThermoInterface
     */
    template <class Base, class T>
    class PyThermoInterface : public Base
    {
        NB_TRAMPOLINE(Base, 1);

    public:
        THERMO_INTERFACE_BASE_TRAMPOLINE();
        THERMO_INTERFACE_TRAMPOLINE(T);
    };

    template <class Base, class T>
    class PyThermoExtendedInterface : public Base
    {
        NB_TRAMPOLINE(Base, 1);

    public:
        THERMO_INTERFACE_BASE_TRAMPOLINE();
        THERMO_EXT_INTERFACE_TRAMPOLINE(T);
    };

#undef THERMO_INTERFACE_BASE_TRAMPOLINE
#undef THERMO_INTERFACE_TRAMPOLINE
#undef THERMO_EXT_INTERFACE_TRAMPOLINE

    void bind_thermo_base(nb::module_& m)
    {
        using namespace nb::literals;

        //  Binding of the abstract/interface class PyThermoInterface
        auto g = nb::class_<PyThermoExtendedInterface<ThermoExtendedInterface<double>, double>>(m, "Thermo");
        g.def(nb::init<>());
        bind_thermo_extended_interface<ThermoExtendedInterface<double>>(g);


        nb::exception<convergence_error>(m, "ConvergenceError", PyExc_RuntimeError);
        nb::exception<domain_error>(m, "DomainError", PyExc_RuntimeError);
    }
}