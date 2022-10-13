// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"

#include "libthermo/real_gas.hpp"


namespace py = pybind11;
using namespace thermo;

class PyRealGas
    : public ThermoInterface<array_t>
    , private RealGas
{
public:
    PyRealGas(double r_)
        : RealGas(r_)
    {
    }

    double Gamma(double t) const override
    {
        return RealGas::Gamma(t);
    }
    array_t Gamma(const array_t& t) const override
    {
        return RealGas::Gamma(t);
    }

    double Cp(double t) const override
    {
        return RealGas::Cp(t);
    }
    array_t Cp(const array_t& t) const override
    {
        return RealGas::Cp(t);
    }

    double H(double t) const override
    {
        return RealGas::H(t);
    }
    array_t H(const array_t& t) const override
    {
        return RealGas::H(t);
    }

    double Phi(double t) const override
    {
        return RealGas::Phi(t);
    }
    array_t Phi(const array_t& t) const override
    {
        return RealGas::Phi(t);
    }

    double PR(double t1, double t2, double eff_poly) const override
    {
        return RealGas::PR(t1, t2, eff_poly);
    }
    array_t PR(const array_t& t1, const array_t& t2, const array_t& eff_poly) const override
    {
        return RealGas::PR(t1, t2, eff_poly);
    }

    double EffPoly(double p1, double t1, double p2, double t2) const override
    {
        return RealGas::EffPoly(p1, t1, p2, t2);
    }
    array_t EffPoly(const array_t& p1,
                    const array_t& t1,
                    const array_t& p2,
                    const array_t& t2) const override
    {
        return RealGas::EffPoly(p1, t1, p2, t2);
    }

    double R() const override
    {
        return RealGas::R();
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
    using array_t = xt::pyarray<double, xt::layout_type::row_major>;

    py::class_<PyRealGas, ThermoInterface<array_t>, std::shared_ptr<PyRealGas>>(m, "RealGas")
        .def(py::init<double>())
        .def_property_readonly("constant", &PyRealGas::R, "Gas constant")
        .def_property_readonly("r", &PyRealGas::R, "Gas constant")

        .def("enthalpy",
             py::overload_cast<double>(&PyRealGas::H, py::const_),
             "Enthalpy",
             py::arg("temperature"))
        .def("enthalpy",
             py::overload_cast<const array_t&>(&PyRealGas::H, py::const_),
             "Enthalpy",
             py::arg("temperature"))
        .def("h",
             py::overload_cast<double>(&PyRealGas::H, py::const_),
             "Enthalpy",
             py::arg("temperature"))
        .def("h",
             py::overload_cast<const array_t&>(&PyRealGas::H, py::const_),
             "Enthalpy",
             py::arg("temperature"))

        .def("entropy",
             py::overload_cast<double>(&PyRealGas::Phi, py::const_),
             "Entropy",
             py::arg("temperature"))
        .def("entropy",
             py::overload_cast<const array_t&>(&PyRealGas::Phi, py::const_),
             "Entropy",
             py::arg("temperature"))
        .def("phi",
             py::overload_cast<double>(&PyRealGas::Phi, py::const_),
             "Entropy",
             py::arg("temperature"))
        .def("phi",
             py::overload_cast<const array_t&>(&PyRealGas::Phi, py::const_),
             "Entropy",
             py::arg("temperature"))
        /*
        .def(
            "dphi",
            [](const PyRealGas& self, const array_t& t1, const array_t& t2) -> array_t {
                return self.dPhi(t1, t2);
            },
            "Entropy delta",
            py::arg("initiale temperature"),
            py::arg("final temperature"))
        */
        .def("specific_heat_ratio",
             py::overload_cast<double>(&PyRealGas::Gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))
        .def("specific_heat_ratio",
             py::overload_cast<const array_t&>(&PyRealGas::Gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))
        .def("gamma",
             py::overload_cast<double>(&PyRealGas::Gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))
        .def("gamma",
             py::overload_cast<const array_t&>(&PyRealGas::Gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))

        .def("specific_heat_pressure",
             py::overload_cast<double>(&PyRealGas::Cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))
        .def("specific_heat_pressure",
             py::overload_cast<const array_t&>(&PyRealGas::Cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp",
             py::overload_cast<double>(&PyRealGas::Cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp",
             py::overload_cast<const array_t&>(&PyRealGas::Cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))

        .def("pressure_ratio",
             py::overload_cast<double, double, double>(&PyRealGas::PR, py::const_),
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pressure_ratio",
             py::overload_cast<const array_t&, const array_t&, const array_t&>(&PyRealGas::PR,
                                                                               py::const_),
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             py::overload_cast<double, double, double>(&PyRealGas::PR, py::const_),
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             py::overload_cast<const array_t&, const array_t&, const array_t&>(&PyRealGas::PR,
                                                                               py::const_),
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))

        .def("eff_poly",
             py::overload_cast<double, double, double, double>(&PyRealGas::EffPoly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("eff_poly",
             py::overload_cast<const array_t&, const array_t&, const array_t&, const array_t&>(
                 &PyRealGas::EffPoly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("polytropic_efficiency",
             py::overload_cast<double, double, double, double>(&PyRealGas::EffPoly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("polytropic_efficiency",
             py::overload_cast<const array_t&, const array_t&, const array_t&, const array_t&>(
                 &PyRealGas::EffPoly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))

        .def("t_from_h", &PyRealGas::TFromH, "Temperature from enthalpy", py::arg("enthalpy"))
        .def("t_from_phi", &PyRealGas::TFromPhi, "Temperature from entropy", py::arg("entropy"))
        .def("t_from_pr",
             &PyRealGas::TFromPR,
             "Temperature from pressure ratio, initial temperature and polytropic efficiency",
             py::arg("pr"),
             py::arg("t1"),
             py::arg("eff_poly"));
}
