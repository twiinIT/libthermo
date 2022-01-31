// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "gas.hpp"

#include "libthermo/ideal_gas.h"
#include "libthermo/real_gas.h"


namespace py = pybind11;
using namespace thermo;

double
PyGas::Gamma(double t) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, gamma, t);
}

array_type
PyGas::Gamma(const array_type& t) const override
{
    PYBIND11_OVERLOAD_PURE(array_type, ThermoInterface, gamma, t);
}

double
PyGas::Cp(double t) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, cp, t);
}

array_type
PyGas::Cp(const array_type& t) const override
{
    PYBIND11_OVERLOAD_PURE(array_type, ThermoInterface, cp, t);
}

double
PyGas::H(double t) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, h, t);
}

array_type
PyGas::H(const array_type& t) const override
{
    PYBIND11_OVERLOAD_PURE(array_type, ThermoInterface, h, t);
}

double
PyGas::Phi(double t) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, phi, t);
}

array_type
PyGas::Phi(const array_type& t) const override
{
    PYBIND11_OVERLOAD_PURE(array_type, ThermoInterface, phi, t);
}

double
PyGas::PR(double t1, double t2, double eff_poly) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, pr, t1, t2, eff_poly);
}

array_type
PyGas::PR(const array_type& t1, const array_type& t2, const array_type& eff_poly) const override
{
    PYBIND11_OVERLOAD_PURE(array_type, ThermoInterface, pr, t1, t2, eff_poly);
}

double
PyGas::EffPoly(double p1, double t1, double p2, double t2) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, eff_poly, p1, t1, p2, t2);
}

array_type
PyGas::EffPoly(const array_type& p1,
               const array_type& t1,
               const array_type& p2,
               const array_type& t2) const override
{
    PYBIND11_OVERLOAD_PURE(array_type, ThermoInterface, eff_poly, p1, t1, p2, t2);
}

double
PyGas::R() const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, r, );
}

double
PyGas::TFromPR(double pr_, double t1_, double eff_poly_) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_pr, pr_, t1_, eff_poly_);
}

double
PyGas::TFromH(double h_) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_h, h_);
}

double
PyGas::TFromPhi(double phi_) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_phi, phi_);
}

double
PyGas::StaticT(double tt_, double mach_) const override
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, static_t, tt_, mach_);
}


void
gas(py::module_& m)
{
    // ==== Binding of the abstract Gas class ==== //
    auto g = py::class_<ThermoInterface<array_type>,
                        PyGas,
                        std::shared_ptr<ThermoInterface<array_type>>>(m, "Gas");
    g.def(py::init<>());
    g.def_property_readonly("constant", &ThermoInterface<array_type>::R, "Gas constant");
    g.def_property_readonly("r", &ThermoInterface<array_type>::R, "Gas constant");
    g.def("enthalpy", &ThermoInterface<array_type>::H, "Enthalpy", py::arg("temperature"));
    g.def("enthalpy", &ThermoInterface<array_type>::H, "Enthalpy", py::arg("temperature"));
    g.def("h", &ThermoInterface<array_type>::H, "Enthalpy", py::arg("temperature"));
    //.def("h", &ThermoInterface<array_type>::H<array_type>, "Enthalpy", py::arg("temperature"))
    g.def("entropy", &ThermoInterface<array_type>::Phi, "Entropy", py::arg("temperature"));
    g.def("phi", &ThermoInterface<array_type>::Phi, "Entropy", py::arg("temperature"));
    //.def("specific_heat_ratio",
    //    &ThermoInterface::SpecificHeatRatio,
    //       "Specific heat ratio",
    //     py::arg("temperature"))
    g.def("gamma",
          py::overload_cast<double>(&ThermoInterface<array_type>::Gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));
    g.def("specific_heat_pressure",
          &ThermoInterface<array_type>::Cp,
          "Specific heat pressure",
          py::arg("temperature"));
    g.def("cp", &ThermoInterface<array_type>::Cp, "Specific heat pressure", py::arg("temperature"));
    g.def("pressure_ratio",
          &ThermoInterface<array_type>::PR,
          "Pressure Ratio",
          py::arg("initiale temperature"),
          py::arg("finale temperature"),
          py::arg("polytropic efficiency"));
    g.def("pr",
          &ThermoInterface<array_type>::PR,
          "Pressure Ratio",
          py::arg("initiale temperature"),
          py::arg("finale temperature"),
          py::arg("polytropic efficiency"));
    g.def("eff_poly",
          &ThermoInterface<array_type>::EffPoly,
          "Polytropic efficiency",
          py::arg("initial pressure"),
          py::arg("initial temperature"),
          py::arg("final pressure"),
          py::arg("final temperature"));
    g.def("polytropic_efficiency",
          &ThermoInterface<array_type>::EffPoly,
          "Polytropic efficiency",
          py::arg("initial pressure"),
          py::arg("initial temperature"),
          py::arg("final pressure"),
          py::arg("final temperature"));
}
