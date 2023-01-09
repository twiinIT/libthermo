// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include "libthermo/exceptions.hpp"
#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"


namespace py = pybind11;
using namespace thermo;

namespace pythermo
{
#define THERMO_INTERFACE_BASE_TRAMPOLINE(Base)                                                                         \
    double r() const override                                                                                          \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(double, Base, r, );                                                                     \
    }

#define THERMO_INTERFACE_TRAMPOLINE(Base, T)                                                                           \
    T gamma(const T& t) const override                                                                                 \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, gamma, t);                                                                     \
    }                                                                                                                  \
    T cp(const T& t) const override                                                                                    \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, cp, t);                                                                        \
    }                                                                                                                  \
    T h(const T& t) const override                                                                                     \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, h, t);                                                                         \
    }                                                                                                                  \
    T phi(const T& t) const override                                                                                   \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, phi, t);                                                                       \
    }                                                                                                                  \
    T pr(const T& t1, const T& t2, const T& eff_poly) const override                                                   \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, pr, t1, t2, eff_poly);                                                         \
    }                                                                                                                  \
    T eff_poly(const T& p1, const T& t1, const T& p2, const T& t2) const override                                      \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, eff_poly, p1, t1, p2, t2);                                                     \
    }

#define THERMO_EXT_INTERFACE_TRAMPOLINE(Base, T)                                                                       \
    T t_f_h(const T& h_, double tol, std::size_t max_iter = 30) const override                                         \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, t_f_h, h_, tol, max_iter);                                                     \
    }                                                                                                                  \
    T t_f_phi(const T& phi_, double tol, std::size_t max_iter = 30) const override                                     \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, t_f_phi, phi_, tol, max_iter);                                                 \
    }                                                                                                                  \
    T t_f_pr(const T& pr_, const T& t1_, const T& eff_poly_, double tol, std::size_t max_iter = 30) const override     \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, t_f_pr, pr_, t1_, eff_poly_, tol, max_iter);                                   \
    }                                                                                                                  \
    T static_t(const T& tt_, const T& mach_, double tol, std::size_t max_iter = 30) const override                     \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, static_t, tt_, mach_, tol, max_iter);                                          \
    }                                                                                                                  \
    T mach_f_wqa(const T& pt, const T& tt, const T& wqa, double tol, std::size_t max_iter = 30) const override         \
    {                                                                                                                  \
        PYBIND11_OVERLOAD_PURE(T, Base, mach_f_wqa, pt, tt, wqa, tol, max_iter);                                       \
    }                                                                                                                  \
    THERMO_INTERFACE_TRAMPOLINE(Base, T)

    /* Defines the C++ trampolines to allow Python users to
     * define their own thermo and use it in code accepting
     * a PyThermoInterface
     */
    template <class Base, class T>
    class PyThermoInterface : public Base
    {
    public:
        THERMO_INTERFACE_BASE_TRAMPOLINE(Base);
        THERMO_INTERFACE_TRAMPOLINE(Base, T);
    };

    template <class Base, class T>
    class PyThermoExtendedInterface : public Base
    {
    public:
        THERMO_INTERFACE_BASE_TRAMPOLINE(Base);
        THERMO_EXT_INTERFACE_TRAMPOLINE(Base, T);
    };

    class PyThermoTrampoline : public PyThermo<array_t>
    {
    public:
        THERMO_INTERFACE_BASE_TRAMPOLINE(PyThermo<array_t>);
        THERMO_INTERFACE_TRAMPOLINE(PyThermo<array_t>, array_t);
        THERMO_EXT_INTERFACE_TRAMPOLINE(PyThermo<array_t>, double);
    };

#undef THERMO_INTERFACE_BASE_TRAMPOLINE
#undef THERMO_INTERFACE_TRAMPOLINE
#undef THERMO_EXT_INTERFACE_TRAMPOLINE

    void thermo_base(py::module_& m)
    {
        using namespace py::literals;

        py::class_<thermo::ThermoExtendedInterface<double>,
                   PyThermoExtendedInterface<ThermoExtendedInterface<double>, double>,
                   std::shared_ptr<thermo::ThermoExtendedInterface<double>>>(m, "ThermoExtendedInterfaceDouble")
            .def(py::init<>());

        py::class_<thermo::ThermoInterface<array_t>,
                   PyThermoInterface<ThermoInterface<array_t>, array_t>,
                   std::shared_ptr<thermo::ThermoInterface<array_t>>>(m, "ThermoInterfaceArr")
            .def(py::init<>());

        //  Binding of the abstract/interface class PyThermoInterface
        auto g = py::class_<PyThermo<array_t>,
                            thermo::ThermoExtendedInterface<double>,
                            thermo::ThermoInterface<array_t>,
                            PyThermoTrampoline,
                            std::shared_ptr<PyThermo<array_t>>>(m, "Thermo");
        g.def(py::init<>());

        g.def_property_readonly("constant", &ThermoInterface<double>::r, "Gas constant");
        g.def("r", &ThermoInterface<double>::r, "Gas constant");

        g.def("enthalpy", &ThermoInterface<double>::h, "Enthalpy", "temperature"_a);
        g.def("enthalpy", &ThermoInterface<array_t>::h, "Enthalpy", "temperature"_a);
        g.def("h", &ThermoInterface<double>::h, "Enthalpy", "temperature"_a);
        g.def("h", &ThermoInterface<array_t>::h, "Enthalpy", "temperature"_a);

        g.def("entropy", &ThermoInterface<double>::phi, "Entropy", "temperature"_a);
        g.def("entropy", &ThermoInterface<array_t>::phi, "Entropy", "temperature"_a);
        g.def("phi", &ThermoInterface<double>::phi, "Entropy", "temperature"_a);
        g.def("phi", &ThermoInterface<array_t>::phi, "Entropy", "temperature"_a);

        g.def("specific_heat_ratio", &ThermoInterface<double>::gamma, "Specific heat ratio", "temperature"_a);
        g.def("specific_heat_ratio", &ThermoInterface<array_t>::gamma, "Specific heat ratio", "temperature"_a);
        g.def("gamma", &ThermoInterface<double>::gamma, "Specific heat ratio", "temperature"_a);
        g.def("gamma", &ThermoInterface<array_t>::gamma, "Specific heat ratio", "temperature"_a);

        g.def("specific_heat_pressure", &ThermoInterface<double>::cp, "Specific heat pressure", "temperature"_a);
        g.def("specific_heat_pressure", &ThermoInterface<array_t>::cp, "Specific heat pressure", "temperature"_a);
        g.def("cp", &ThermoInterface<double>::cp, "Specific heat pressure", "temperature"_a);
        g.def("cp", &ThermoInterface<array_t>::cp, "Specific heat pressure", "temperature"_a);

        g.def("pressure_ratio", &ThermoInterface<double>::pr, "Pressure Ratio", "t1"_a, "t2"_a, "eff_poly"_a);
        g.def("pressure_ratio", &ThermoInterface<array_t>::pr, "Pressure Ratio", "t1"_a, "t2"_a, "eff_poly"_a);
        g.def("pr", &ThermoInterface<double>::pr, "Pressure Ratio", "t1"_a, "t2"_a, "eff_poly"_a);
        g.def("pr", &ThermoInterface<array_t>::pr, "Pressure Ratio", "t1"_a, "t2"_a, "eff_poly"_a);

        g.def("polytropic_efficiency",
              &ThermoInterface<double>::eff_poly,
              "Polytropic efficiency",
              "p1"_a,
              "t1"_a,
              "p2"_a,
              "t2"_a);
        g.def("polytropic_efficiency",
              &ThermoInterface<array_t>::eff_poly,
              "Polytropic efficiency",
              "p1"_a,
              "t1"_a,
              "p2"_a,
              "t2"_a);
        g.def("eff_poly", &ThermoInterface<double>::eff_poly, "Polytropic efficiency", "p1"_a, "t1"_a, "p2"_a, "t2"_a);
        g.def("eff_poly", &ThermoInterface<array_t>::eff_poly, "Polytropic efficiency", "p1"_a, "t1"_a, "p2"_a, "t2"_a);

        g.def("static_t",
              &ThermoExtendedInterface<double>::static_t,
              "Static temperature",
              "tt"_a,
              "mach"_a,
              "tol"_a,
              "max_iter"_a = 30);

        g.def("t_f_h",
              &ThermoExtendedInterface<double>::t_f_h,
              "Temperature from enthalpy",
              "h"_a,
              "tol"_a,
              "max_iter"_a = 30);

        g.def("t_f_phi",
              &ThermoExtendedInterface<double>::t_f_phi,
              "Temperature from phi function",
              "h"_a,
              "tol"_a,
              "max_iter"_a = 30);

        g.def("t_f_pr",
              &ThermoExtendedInterface<double>::t_f_pr,
              "Temperature from pressure ratio",
              "pr"_a,
              "t1"_a,
              "eff_poly"_a,
              "tol"_a,
              "max_iter"_a = 30);

        g.def("mach_f_wqa",
              &ThermoExtendedInterface<double>::mach_f_wqa,
              "Mach number from specific mass flow",
              "pt"_a,
              "tt"_a,
              "wqa"_a,
              "tol"_a,
              "max_iter"_a = 30);

        py::register_exception<convergence_error>(m, "ConvergenceError", PyExc_RuntimeError);
        py::register_exception<domain_error>(m, "DomainError", PyExc_RuntimeError);
    }
}