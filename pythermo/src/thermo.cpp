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

double
PyThermo::gamma(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, gamma, t);
}

array_t
PyThermo::gamma(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_gamma, t);
}

double
PyThermo::cp(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, cp, t);
}

array_t
PyThermo::cp(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_cp, t);
}

double
PyThermo::h(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, h, t);
}

array_t
PyThermo::h(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_h, t);
}

double
PyThermo::phi(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, phi, t);
}

array_t
PyThermo::phi(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_phi, t);
}

double
PyThermo::pr(double t1, double t2, double eff_poly) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, pr, t1, t2, eff_poly);
}

array_t
PyThermo::pr(const array_t& t1, const array_t& t2, const array_t& eff_poly) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_pr, t1, t2, eff_poly);
}

double
PyThermo::eff_poly(double p1, double t1, double p2, double t2) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, eff_poly, p1, t1, p2, t2);
}

array_t
PyThermo::eff_poly(const array_t& p1, const array_t& t1, const array_t& p2, const array_t& t2) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_eff_poly, p1, t1, p2, t2);
}

double
PyThermo::r() const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, r, );
}

double
PyThermo::t_f_pr(double pr_, double t1_, double eff_poly_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_pr, pr_, t1_, eff_poly_, tol, max_iter);
}

double
PyThermo::t_f_h(double h_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_h, h_, tol, max_iter);
}

double
PyThermo::t_f_phi(double phi_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_phi, phi_, tol, max_iter);
}

double
PyThermo::static_t(double tt_, double mach_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, static_t, tt_, mach_, tol, max_iter);
}

double
PyThermo::mach_f_wqa(double pt_, double tt_, double wqa_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, mach_f_wqa, pt_, tt_, wqa_, tol, max_iter);
}

void
thermo_base(py::module_& m)
{
    // ==== Binding of the abstract Gas class ==== //
    auto g
        = py::class_<ThermoInterface<array_t>, PyThermo, std::shared_ptr<ThermoInterface<array_t>>>(
            m, "Thermo");
    g.def(py::init<>());

    g.def_property_readonly("constant", &ThermoInterface<array_t>::r, "Gas constant");
    g.def_property_readonly("r", &ThermoInterface<array_t>::r, "Gas constant");

    g.def("enthalpy",
          py::overload_cast<double>(&ThermoInterface<array_t>::h, py::const_),
          "Enthalpy",
          py::arg("temperature"));
    g.def("enthalpy",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::h, py::const_),
          "Enthalpy",
          py::arg("temperature"));
    g.def("h",
          py::overload_cast<double>(&ThermoInterface<array_t>::h, py::const_),
          "Enthalpy",
          py::arg("temperature"));
    g.def("h",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::h, py::const_),
          "Enthalpy",
          py::arg("temperature"));

    g.def("entropy",
          py::overload_cast<double>(&ThermoInterface<array_t>::phi, py::const_),
          "Entropy",
          py::arg("temperature"));
    g.def("phi",
          py::overload_cast<double>(&ThermoInterface<array_t>::phi, py::const_),
          "Entropy",
          py::arg("temperature"));

    //.def("specific_heat_ratio",
    //    &ThermoInterface::specific_heat_ratio,
    //       "Specific heat ratio",
    //     py::arg("temperature"))
    g.def("specific_heat_ratio",
          py::overload_cast<double>(&ThermoInterface<array_t>::gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));
    g.def("specific_heat_ratio",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));
    g.def("gamma",
          py::overload_cast<double>(&ThermoInterface<array_t>::gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));
    g.def("gamma",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));

    g.def("specific_heat_pressure",
          py::overload_cast<double>(&ThermoInterface<array_t>::cp, py::const_),
          "Specific heat pressure",
          py::arg("temperature"));
    g.def("cp",
          py::overload_cast<double>(&ThermoInterface<array_t>::cp, py::const_),
          "Specific heat pressure",
          py::arg("temperature"));

    g.def("pressure_ratio",
          py::overload_cast<double, double, double>(&ThermoInterface<array_t>::pr, py::const_),
          "Pressure Ratio",
          py::arg("initiale temperature"),
          py::arg("finale temperature"),
          py::arg("polytropic efficiency"));
    g.def("pr",
          py::overload_cast<double, double, double>(&ThermoInterface<array_t>::pr, py::const_),
          "Pressure Ratio",
          py::arg("initiale temperature"),
          py::arg("finale temperature"),
          py::arg("polytropic efficiency"));

    g.def("eff_poly",
          py::overload_cast<double, double, double, double>(&ThermoInterface<array_t>::eff_poly,
                                                            py::const_),
          "Polytropic efficiency",
          py::arg("initial pressure"),
          py::arg("initial temperature"),
          py::arg("final pressure"),
          py::arg("final temperature"));
    g.def("polytropic_efficiency",
          py::overload_cast<double, double, double, double>(&ThermoInterface<array_t>::eff_poly,
                                                            py::const_),
          "Polytropic efficiency",
          py::arg("initial pressure"),
          py::arg("initial temperature"),
          py::arg("final pressure"),
          py::arg("final temperature"));

    py::register_exception<convergence_error>(m, "ConvergenceError", PyExc_RuntimeError);
    py::register_exception<domain_error>(m, "DomainError", PyExc_RuntimeError);
}
