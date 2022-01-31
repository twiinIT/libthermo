// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "gas.hpp"

#include "libthermo/real_gas.hpp"

#include <pybind11/pybind11.h>
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <iostream>
#include <chrono>
#include <functional>


class PyRealGas
    : public ThermoInterface<array_type>
    , public RealGas
{
public:
    PyRealGas(double r_)
        : RealGas(r_)
    {
    }
    using RealGas::Cp;
    using RealGas::EffPoly;
    using RealGas::Gamma;
    using RealGas::H;
    using RealGas::Phi;
    using RealGas::PR;

    double Gamma(double t) const override
    {
        return RealGas::Gamma(t);
    }
    array_type Gamma(const array_type& t) const override
    {
        return RealGas::Gamma(t);
    }
    double Cp(double t) const override
    {
        return RealGas::Cp(t);
    }
    double H(double t) const override
    {
        return RealGas::H(t);
    }
    double Phi(double t) const override
    {
        return RealGas::Phi(t);
    }
    double R() const override
    {
        return RealGas::R();
    }
    double PR(double t1, double t2, double eff_poly) const override
    {
        return RealGas::PR(t1, t2, eff_poly);
    }
    double EffPoly(double p1, double t1, double p2, double t2) const override
    {
        return RealGas::EffPoly(p1, t1, p2, t2);
    }
    double StaticT(double t, double mach) const override
    {
        return RealGas::StaticT(t, mach);
    }
    double TFromPR(double pr, double t1, double eff_poly) const override
    {
        return RealGas::TFromPR(pr, t1, eff_poly);
    }
    double TFromH(double h) const override
    {
        return RealGas::TFromH(h);
    }
    double TFromPhi(double phi) const override
    {
        return RealGas::TFromPhi(phi);
    }
};


void
real_gas(py::module_& m)
{
    using array_type = xt::pyarray<double, xt::layout_type::row_major>;

    py::class_<PyRealGas, ThermoInterface<array_type>, std::shared_ptr<PyRealGas>>(m, "RealGas")
        .def(py::init<double>())
        .def_property_readonly("constant", &RealGas::Constant, "Gas constant")
        .def_property_readonly("r", &RealGas::R, "Gas constant")
        .def("enthalpy", &RealGas::Enthalpy<double>, "Enthalpy", py::arg("temperature"))
        .def("h", &RealGas::H<double>, "Enthalpy", py::arg("temperature"))
        .def(
            "h",
            [](const PyRealGas& self, const array_type& t1) -> array_type { return self.H(t1); },
            "Enthalpy",
            py::arg("temperature"))
        .def("entropy", &RealGas::Entropy<double>, "Entropy", py::arg("temperature"))
        .def("phi", &RealGas::Phi<double>, "Entropy", py::arg("temperature"))
        .def(
            "phi",
            [](const PyRealGas& self, const array_type& t1) -> array_type { return self.Phi(t1); },
            "Entropy",
            py::arg("temperature"))
        .def(
            "dphi",
            [](const PyRealGas& self, const array_type& t1, const array_type& t2) -> array_type {
                return self.dPhi(t1, t2);
            },
            "Entropy delta",
            py::arg("initiale temperature"),
            py::arg("final temperature"))
        .def("specific_heat_ratio",
             &RealGas::SpecificHeatRatio<double>,
             "Specific heat ratio",
             py::arg("temperature"))
        .def(
            "gamma",
            [](const PyRealGas& self, const array_type& t1) -> array_type {
                return self.Gamma(t1);
            },
            "Specific heat ratio",
            py::arg("temperature"))
        .def("gamma", &RealGas::Gamma<double>, "Specific heat ratio", py::arg("temperature"))
        .def("specific_heat_pressure",
             &RealGas::SpecificHeatPressure<double>,
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp", &RealGas::Cp<double>, "Specific heat pressure", py::arg("temperature"))
        .def(
            "cp",
            [](const PyRealGas& self, const array_type& t) -> array_type { return self.Cp(t); },
            "Specific heat pressure",
            py::arg("temperature"))
        .def("pressure_ratio",
             &RealGas::PressureRatio<double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             &RealGas::PR<double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def(
            "pr",
            [](const PyRealGas& self,
               const array_type& t1,
               const array_type& t2,
               const array_type& eff_poly) -> array_type { return self.PR(t1, t2, eff_poly); },
            "Pressure Ratio",
            py::arg("initiale temperature"),
            py::arg("finale temperature"),
            py::arg("polytropic efficiency"))
        .def("pr",
             &RealGas::PR<array_type, double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("eff_poly",
             &RealGas::EffPoly<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def(
            "eff_poly",
            [](const PyRealGas& self,
               const array_type& p1,
               const array_type& t1,
               const array_type& p2,
               const array_type& t2) -> array_type { return self.EffPoly(p1, t1, p2, t2); },
            "Polytropic efficiency",
            py::arg("initial pressure"),
            py::arg("initial temperature"),
            py::arg("final pressure"),
            py::arg("final temperature"))
        .def("polytropic_efficiency",
             &RealGas::PolytropicEfficiency<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("t_from_h", &RealGas::TFromH<double>, "Temperature from enthalpy", py::arg("enthalpy"))
        .def("t_from_phi",
             &RealGas::TFromPhi<double>,
             "Temperature from entropy",
             py::arg("entropy"))
        .def("t_from_pr",
             &RealGas::TFromPR<double>,
             "Temperature from pressure ratio, initial temperature and polytropic efficiency",
             py::arg("pr"),
             py::arg("t1"),
             py::arg("eff_poly"));
}
