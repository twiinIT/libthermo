// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "libthermo/ideal_gas.h"
#include "libthermo/real_gas.h"

#include <pybind11/pybind11.h>
#define FORCE_IMPORT_ARRAY
#include "xtensor-python/pytensor.hpp"
#include "xtensor-python/pyarray.hpp"

#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <iostream>
#include <chrono>
#include <functional>


namespace py = pybind11;
using namespace thermo;
using array_type = xt::pyarray<double, xt::layout_type::row_major>;

class PyGas : public GasInterface<array_type>
{
public:
    using GasInterface::GasInterface;

    double Gamma(double t) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, gamma, t);
    }
    array_type Gamma(const array_type& t) const override
    {
        PYBIND11_OVERLOAD_PURE(array_type, GasInterface, gamma, t);
    }
    double Cp(double t) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, cp, t);
    }
    double H(double t) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, h, t);
    }
    double Phi(double t) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, phi, t);
    }
    double R() const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, r, );
    }
    double PR(double t1, double t2, double eff_poly) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, pr, t1, t2, eff_poly);
    }
    double EffPoly(double p1, double t1, double p2, double t2) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, eff_poly, p1, t1, p2, t2);
    }

    double TFromPR(double pr_, double t1_, double eff_poly_) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, T_from_pr, pr_, t1_, eff_poly_);
    }

    double TFromH(double h_) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, T_from_h, h_);
    }

    double TFromPhi(double phi_) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, T_from_phi, phi_);
    }

    double StaticT(double tt_, double mach_) const override
    {
        PYBIND11_OVERLOAD_PURE(double, GasInterface, static_T, tt_, mach_);
    }
};


class PyIdealGas
    : public GasInterface<array_type>
    , public IdealGas
{
public:
    PyIdealGas(double r_, double cp_)
        : IdealGas(r_, cp_)
    {
    }
    using IdealGas::Cp;
    using IdealGas::EffPoly;
    using IdealGas::Gamma;
    using IdealGas::H;
    using IdealGas::Phi;
    using IdealGas::PR;

    double Gamma(double t) const override
    {
        return IdealGas::Gamma(t);
    }
    array_type Gamma(const array_type& t) const override
    {
        return IdealGas::Gamma(t);
    }
    double Cp(double t) const override
    {
        return IdealGas::Cp(t);
    }
    double H(double t) const override
    {
        return IdealGas::H(t);
    }
    double Phi(double t) const override
    {
        return IdealGas::Phi(t);
    }
    double R() const override
    {
        return IdealGas::R();
    }
    double PR(double t1, double t2, double eff_poly) const override
    {
        return IdealGas::PR(t1, t2, eff_poly);
    }
    double EffPoly(double p1, double t1, double p2, double t2) const override
    {
        return IdealGas::EffPoly(p1, t1, p2, t2);
    }
    double StaticT(double t, double mach) const override
    {
        return IdealGas::StaticT(t, mach);
    }
    double TFromPR(double pr, double t1, double eff_poly) const override
    {
        return IdealGas::TFromPR(pr, t1, eff_poly);
    }
    double TFromH(double h) const override
    {
        return IdealGas::TFromH(h);
    }
    double TFromPhi(double phi) const override
    {
        return IdealGas::TFromPhi(phi);
    }
};

class PyRealGas
    : public GasInterface<array_type>
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
print_h(GasInterface<array_type>* gas)
{
    double h = gas->H(288.15);
    std::cout << "h: " << h << std::endl;
}


PYBIND11_MODULE(pythermo, m)
{
    m.doc() = "Python bindings of the C++ implementation of the thermodynamic library 'libthermo'.";

    xt::import_numpy();

    using array_type = xt::pyarray<double, xt::layout_type::row_major>;

    // ==== Binding of the abstract Gas class ==== //
    py::class_<GasInterface<array_type>, PyGas, std::shared_ptr<GasInterface<array_type>>>(m, "Gas")
        .def(py::init<>())
        //.def_property_readonly("constant", &GasInterface::Constant, "Gas constant")
        .def_property_readonly("r", &GasInterface<array_type>::R, "Gas constant")
        //.def("enthalpy", &GasInterface<array_type>::Enthalpy, "Enthalpy", py::arg("temperature"))
        //.def("h", &GasInterface<array_type>::H<double>, "Enthalpy", py::arg("temperature"))
        //.def("h", &GasInterface<array_type>::H<array_type>, "Enthalpy", py::arg("temperature"))
        //.def("entropy", &GasInterface::Entropy, "Entropy", py::arg("temperature"))
        .def("phi", &GasInterface<array_type>::Phi, "Entropy", py::arg("temperature"))
        //.def("specific_heat_ratio",
        //    &GasInterface::SpecificHeatRatio,
        //       "Specific heat ratio",
        //     py::arg("temperature"))
        .def("gamma",
             py::overload_cast<double>(&GasInterface<array_type>::Gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))
        //.def("specific_heat_pressure",
        //     &GasInterface::SpecificHeatPressure,
        //     "Specific heat pressure",
        //     py::arg("temperature"))
        .def("cp", &GasInterface<array_type>::Cp, "Specific heat pressure", py::arg("temperature"))
        .def("pressure_ratio",
             &GasInterface<array_type>::PR,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             &GasInterface<array_type>::PR,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("eff_poly",
             &GasInterface<array_type>::EffPoly,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        /*.def("polytropic_efficiency",
             &GasInterface::PolytropicEfficiency,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))*/
        ;

    m.def("print_h", &print_h);

    py::class_<PyIdealGas, GasInterface<array_type>, std::shared_ptr<PyIdealGas>>(m, "IdealGas")
        .def(py::init<double, double>())
        .def_property_readonly("constant", &IdealGas::Constant, "Gas constant")
        .def_property_readonly("r", &IdealGas::R, "Gas constant")
        .def("enthalpy", &IdealGas::Enthalpy<double>, "Enthalpy", py::arg("temperature"))
        .def("h", &IdealGas::H<double>, "Enthalpy", py::arg("temperature"))
        .def("entropy", &IdealGas::Entropy<double>, "Entropy", py::arg("temperature"))
        .def("phi", &IdealGas::Phi<double>, "Entropy", py::arg("temperature"))
        .def("phi", &IdealGas::Phi<array_type>, "Entropy", py::arg("temperature"))
        .def("specific_heat_ratio",
             &IdealGas::SpecificHeatRatio<double>,
             "Specific heat ratio",
             py::arg("temperature"))
        .def("gamma", &IdealGas::Gamma<double>, "Specific heat ratio", py::arg("temperature"))
        .def("gamma", &IdealGas::Gamma<array_type>, "Specific heat ratio", py::arg("temperature"))
        .def("specific_heat_pressure",
             &IdealGas::SpecificHeatPressure<double>,
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp", &IdealGas::Cp<double>, "Specific heat pressure", py::arg("temperature"))
        .def("cp", &IdealGas::Cp<array_type>, "Specific heat pressure", py::arg("temperature"))
        .def("pressure_ratio",
             &IdealGas::PressureRatio<double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             &IdealGas::PR<double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             &IdealGas::PR<array_type, array_type>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             &IdealGas::PR<array_type, double>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("eff_poly",
             &IdealGas::EffPoly<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("eff_poly",
             &IdealGas::PolytropicEfficiency<array_type>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("polytropic_efficiency",
             &IdealGas::PolytropicEfficiency<double>,
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def(
            "t_from_h", &IdealGas::TFromH<double>, "Temperature from enthalpy", py::arg("enthalpy"))
        .def("t_from_phi",
             &IdealGas::TFromPhi<double>,
             "Temperature from entropy",
             py::arg("entropy"))
        .def("t_from_pr",
             &IdealGas::TFromPR<double>,
             "Temperature from pressure ratio, initial temperature and polytropic efficiency",
             py::arg("pr"),
             py::arg("t1"),
             py::arg("eff_poly"));

    py::class_<PyRealGas, GasInterface<array_type>, std::shared_ptr<PyRealGas>>(m, "RealGas")
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
