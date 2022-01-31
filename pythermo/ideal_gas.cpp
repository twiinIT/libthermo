// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "gas.hpp"

#include "libthermo/thermo.hpp"
#include "libthermo/ideal_gas.hpp"


namespace py = pybind11;
using namespace thermo;


class PyIdealGas
    : public ThermoInterface<xt::pyarray<double, xt::layout_type::row_major>>
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
    array_t Gamma(const array_t& t) const override
    {
        return IdealGas::Gamma(t);
    }

    double Cp(double t) const override
    {
        return IdealGas::Cp(t);
    }
    array_t Cp(const array_t& t) const override
    {
        return IdealGas::Cp(t);
    }

    double H(double t) const override
    {
        return IdealGas::H(t);
    }
    array_t H(const array_t& t) const override
    {
        return IdealGas::H(t);
    }

    double Phi(double t) const override
    {
        return IdealGas::Phi(t);
    }
    array_t Phi(const array_t& t) const override
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
    array_t PR(const array_t& t1, const array_t& t2, const array_t& eff_poly) const override
    {
        return IdealGas::PR(t1, t2, eff_poly);
    }

    double EffPoly(double p1, double t1, double p2, double t2) const override
    {
        return IdealGas::EffPoly(p1, t1, p2, t2);
    }
    array_t EffPoly(const array_t& p1,
                    const array_t& t1,
                    const array_t& p2,
                    const array_t& t2) const override
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


void
ideal_gas(py::module_& m)
{
    py::class_<PyIdealGas, ThermoInterface<array_t>, std::shared_ptr<PyIdealGas>>(m, "IdealGas")
        .def(py::init<double, double>())
        .def_property_readonly("constant", &IdealGas::Constant, "Gas constant")
        .def_property_readonly("r", &IdealGas::R, "Gas constant")
        .def("enthalpy", &IdealGas::Enthalpy<double>, "Enthalpy", py::arg("temperature"))
        .def("h", &IdealGas::H<double>, "Enthalpy", py::arg("temperature"))
        .def("entropy", &IdealGas::Entropy<double>, "Entropy", py::arg("temperature"))
        .def("phi", &IdealGas::Phi<double>, "Entropy", py::arg("temperature"))
        .def("phi", &IdealGas::Phi<array_t>, "Entropy", py::arg("temperature"))
        .def("specific_heat_ratio",
             &IdealGas::SpecificHeatRatio<double>,
             "Specific heat ratio",
             py::arg("temperature"))
        .def("gamma", &IdealGas::Gamma<double>, "Specific heat ratio", py::arg("temperature"))
        .def("gamma", &IdealGas::Gamma<array_t>, "Specific heat ratio", py::arg("temperature"))
        .def("specific_heat_pressure",
             &IdealGas::SpecificHeatPressure<double>,
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp", &IdealGas::Cp<double>, "Specific heat pressure", py::arg("temperature"))
        .def("cp", &IdealGas::Cp<array_t>, "Specific heat pressure", py::arg("temperature"))
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
             &IdealGas::PR<array_t, array_t>,
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             &IdealGas::PR<array_t, double>,
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
             &IdealGas::PolytropicEfficiency<array_t>,
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
}
