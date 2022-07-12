// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include "libthermo/ideal_gas.hpp"
#include "libthermo/real_gas.hpp"


namespace py = pybind11;
using namespace thermo;

double
PyThermo::Gamma(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, gamma, t);
}

array_t
PyThermo::Gamma(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_gamma, t);
}

double
PyThermo::Cp(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, cp, t);
}

array_t
PyThermo::Cp(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_cp, t);
}

double
PyThermo::H(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, h, t);
}

array_t
PyThermo::H(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_h, t);
}

double
PyThermo::Phi(double t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, phi, t);
}

array_t
PyThermo::Phi(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_phi, t);
}

double
PyThermo::PR(double t1, double t2, double eff_poly) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, pr, t1, t2, eff_poly);
}

array_t
PyThermo::PR(const array_t& t1, const array_t& t2, const array_t& eff_poly) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_pr, t1, t2, eff_poly);
}

double
PyThermo::EffPoly(double p1, double t1, double p2, double t2) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, eff_poly, p1, t1, p2, t2);
}

array_t
PyThermo::EffPoly(const array_t& p1, const array_t& t1, const array_t& p2, const array_t& t2) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface, v_eff_poly, p1, t1, p2, t2);
}

double
PyThermo::R() const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, r, );
}

double
PyThermo::TFromPR(double pr_, double t1_, double eff_poly_) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_pr, pr_, t1_, eff_poly_);
}

double
PyThermo::TFromH(double h_) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_h, h_);
}

double
PyThermo::TFromPhi(double phi_) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, t_from_phi, phi_);
}

double
PyThermo::StaticT(double tt_, double mach_) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface, static_t, tt_, mach_);
}


void
thermo_base(py::module_& m)
{
    // ==== Binding of the abstract Gas class ==== //
    auto g
        = py::class_<ThermoInterface<array_t>, PyThermo, std::shared_ptr<ThermoInterface<array_t>>>(
            m, "Thermo");
    g.def(py::init<>());

    g.def_property_readonly("constant", &ThermoInterface<array_t>::R, "Gas constant");
    g.def_property_readonly("r", &ThermoInterface<array_t>::R, "Gas constant");

    g.def("enthalpy",
          py::overload_cast<double>(&ThermoInterface<array_t>::H, py::const_),
          "Enthalpy",
          py::arg("temperature"));
    g.def("enthalpy",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::H, py::const_),
          "Enthalpy",
          py::arg("temperature"));
    g.def("h",
          py::overload_cast<double>(&ThermoInterface<array_t>::H, py::const_),
          "Enthalpy",
          py::arg("temperature"));
    g.def("h",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::H, py::const_),
          "Enthalpy",
          py::arg("temperature"));

    g.def("entropy",
          py::overload_cast<double>(&ThermoInterface<array_t>::Phi, py::const_),
          "Entropy",
          py::arg("temperature"));
    g.def("phi",
          py::overload_cast<double>(&ThermoInterface<array_t>::Phi, py::const_),
          "Entropy",
          py::arg("temperature"));

    //.def("specific_heat_ratio",
    //    &ThermoInterface::SpecificHeatRatio,
    //       "Specific heat ratio",
    //     py::arg("temperature"))
    g.def("specific_heat_ratio",
          py::overload_cast<double>(&ThermoInterface<array_t>::Gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));
    g.def("specific_heat_ratio",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::Gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));
    g.def("gamma",
          py::overload_cast<double>(&ThermoInterface<array_t>::Gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));
    g.def("gamma",
          py::overload_cast<const array_t&>(&ThermoInterface<array_t>::Gamma, py::const_),
          "Specific heat ratio",
          py::arg("temperature"));

    g.def("specific_heat_pressure",
          py::overload_cast<double>(&ThermoInterface<array_t>::Cp, py::const_),
          "Specific heat pressure",
          py::arg("temperature"));
    g.def("cp",
          py::overload_cast<double>(&ThermoInterface<array_t>::Cp, py::const_),
          "Specific heat pressure",
          py::arg("temperature"));

    g.def("pressure_ratio",
          py::overload_cast<double, double, double>(&ThermoInterface<array_t>::PR, py::const_),
          "Pressure Ratio",
          py::arg("initiale temperature"),
          py::arg("finale temperature"),
          py::arg("polytropic efficiency"));
    g.def("pr",
          py::overload_cast<double, double, double>(&ThermoInterface<array_t>::PR, py::const_),
          "Pressure Ratio",
          py::arg("initiale temperature"),
          py::arg("finale temperature"),
          py::arg("polytropic efficiency"));

    g.def("eff_poly",
          py::overload_cast<double, double, double, double>(&ThermoInterface<array_t>::EffPoly,
                                                            py::const_),
          "Polytropic efficiency",
          py::arg("initial pressure"),
          py::arg("initial temperature"),
          py::arg("final pressure"),
          py::arg("final temperature"));
    g.def("polytropic_efficiency",
          py::overload_cast<double, double, double, double>(&ThermoInterface<array_t>::EffPoly,
                                                            py::const_),
          "Polytropic efficiency",
          py::arg("initial pressure"),
          py::arg("initial temperature"),
          py::arg("final pressure"),
          py::arg("final temperature"));
}
