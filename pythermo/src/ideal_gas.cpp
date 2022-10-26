// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"
#include "libthermo/ideal_gas.hpp"


namespace py = pybind11;
using namespace thermo;

class PyIdealGas
    : public ThermoInterface<array_t>
    , private IdealGas
{
public:
    PyIdealGas(double r_, double cp_)
        : IdealGas(r_, cp_)
    {
    }

    double gamma(double t) const override
    {
        return IdealGas::gamma(t);
    }
    array_t gamma(const array_t& t) const override
    {
        return IdealGas::gamma(t);
    }

    double cp(double t) const override
    {
        return IdealGas::cp(t);
    }
    array_t cp(const array_t& t) const override
    {
        return IdealGas::cp(t);
    }

    double h(double t) const override
    {
        return IdealGas::h(t);
    }
    array_t h(const array_t& t) const override
    {
        return IdealGas::h(t);
    }

    double phi(double t) const override
    {
        return IdealGas::phi(t);
    }
    array_t phi(const array_t& t) const override
    {
        return IdealGas::phi(t);
    }

    double r() const override
    {
        return IdealGas::r();
    }

    double pr(double t1, double t2, double eff_poly) const override
    {
        return IdealGas::pr(t1, t2, eff_poly);
    }
    array_t pr(const array_t& t1, const array_t& t2, const array_t& eff_poly) const override
    {
        return IdealGas::pr(t1, t2, eff_poly);
    }

    double eff_poly(double p1, double t1, double p2, double t2) const override
    {
        return IdealGas::eff_poly(p1, t1, p2, t2);
    }
    array_t eff_poly(const array_t& p1,
                     const array_t& t1,
                     const array_t& p2,
                     const array_t& t2) const override
    {
        return IdealGas::eff_poly(p1, t1, p2, t2);
    }

    double static_t(double t, double mach, double = 1e-8, std::size_t = 1) const override
    {
        return IdealGas::static_t(t, mach);
    }
    double t_f_pr(
        double pr, double t1, double eff_poly, double = 1e-8, std::size_t = 1) const override
    {
        return IdealGas::t_f_pr(pr, t1, eff_poly);
    }
    double t_f_h(double h, double = 1e-8, std::size_t = 1) const override
    {
        return IdealGas::t_f_h(h);
    }
    double t_f_phi(double phi, double = 1e-8, std::size_t = 1) const override
    {
        return IdealGas::t_f_phi(phi);
    }
    double mach_f_wqa(
        double pt, double tt, double wqa, double = 1e-8, std::size_t = 1) const override
    {
        return IdealGas::mach_f_wqa(pt, tt, wqa);
    }
};


void
ideal_gas(py::module_& m)
{
    using array_t = xt::pyarray<double, xt::layout_type::row_major>;

    py::class_<PyIdealGas, ThermoInterface<array_t>, std::shared_ptr<PyIdealGas>>(m, "IdealGas")
        .def(py::init<double, double>())
        .def_property_readonly("constant", &PyIdealGas::r, "Gas constant")
        .def_property_readonly("r", &PyIdealGas::r, "Gas constant")

        .def("enthalpy",
             py::overload_cast<double>(&PyIdealGas::h, py::const_),
             "Enthalpy",
             py::arg("temperature"))
        .def("enthalpy",
             py::overload_cast<const array_t&>(&PyIdealGas::h, py::const_),
             "Enthalpy",
             py::arg("temperature"))
        .def("h",
             py::overload_cast<double>(&PyIdealGas::h, py::const_),
             "Enthalpy",
             py::arg("temperature"))
        .def("h",
             py::overload_cast<const array_t&>(&PyIdealGas::h, py::const_),
             "Enthalpy",
             py::arg("temperature"))

        .def("entropy",
             py::overload_cast<double>(&PyIdealGas::phi, py::const_),
             "Entropy",
             py::arg("temperature"))
        .def("entropy",
             py::overload_cast<const array_t&>(&PyIdealGas::phi, py::const_),
             "Entropy",
             py::arg("temperature"))
        .def("phi",
             py::overload_cast<double>(&PyIdealGas::phi, py::const_),
             "Entropy",
             py::arg("temperature"))
        .def("phi",
             py::overload_cast<const array_t&>(&PyIdealGas::phi, py::const_),
             "Entropy",
             py::arg("temperature"))

        .def("specific_heat_ratio",
             py::overload_cast<double>(&PyIdealGas::gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))
        .def("specific_heat_ratio",
             py::overload_cast<const array_t&>(&PyIdealGas::gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))
        .def("gamma",
             py::overload_cast<double>(&PyIdealGas::gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))
        .def("gamma",
             py::overload_cast<const array_t&>(&PyIdealGas::gamma, py::const_),
             "Specific heat ratio",
             py::arg("temperature"))

        .def("specific_heat_pressure",
             py::overload_cast<double>(&PyIdealGas::cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))
        .def("specific_heat_pressure",
             py::overload_cast<const array_t&>(&PyIdealGas::cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp",
             py::overload_cast<double>(&PyIdealGas::cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))
        .def("cp",
             py::overload_cast<const array_t&>(&PyIdealGas::cp, py::const_),
             "Specific heat pressure",
             py::arg("temperature"))

        .def("pressure_ratio",
             py::overload_cast<double, double, double>(&PyIdealGas::pr, py::const_),
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             py::overload_cast<double, double, double>(&PyIdealGas::pr, py::const_),
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))
        .def("pr",
             py::overload_cast<const array_t&, const array_t&, const array_t&>(&PyIdealGas::pr,
                                                                               py::const_),
             "Pressure Ratio",
             py::arg("initiale temperature"),
             py::arg("finale temperature"),
             py::arg("polytropic efficiency"))

        .def("eff_poly",
             py::overload_cast<double, double, double, double>(&PyIdealGas::eff_poly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("eff_poly",
             py::overload_cast<const array_t&, const array_t&, const array_t&, const array_t&>(
                 &PyIdealGas::eff_poly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("polytropic_efficiency",
             py::overload_cast<double, double, double, double>(&PyIdealGas::eff_poly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))
        .def("polytropic_efficiency",
             py::overload_cast<const array_t&, const array_t&, const array_t&, const array_t&>(
                 &PyIdealGas::eff_poly, py::const_),
             "Polytropic efficiency",
             py::arg("initial pressure"),
             py::arg("initial temperature"),
             py::arg("final pressure"),
             py::arg("final temperature"))

        .def("t_from_h",
             &PyIdealGas::t_f_h,
             "Temperature from enthalpy",
             py::arg("enthalpy"),
             py::arg("tol"),
             py::arg("max_iter"))
        .def("t_from_phi",
             &PyIdealGas::t_f_phi,
             "Temperature from entropy",
             py::arg("entropy"),
             py::arg("tol"),
             py::arg("max_iter"))
        .def("t_from_pr",
             &PyIdealGas::t_f_pr,
             "Temperature from pressure ratio, initial temperature and polytropic efficiency",
             py::arg("pr"),
             py::arg("t1"),
             py::arg("eff_poly"),
             py::arg("tol"),
             py::arg("max_iter"));
}
