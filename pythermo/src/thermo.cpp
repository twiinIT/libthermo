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

/* Defines the C++ trampolines to allow Python users to
 * define their own thermo and use it in code accepting
 * a PyThermoInterface
 */
class PyThermo : public PyThermoInterface
{
public:
    double gamma(const double& t) const override;
    array_t gamma(const array_t& t) const override;

    double cp(const double& t) const override;
    array_t cp(const array_t& t) const;

    double h(const double& t) const override;
    array_t h(const array_t& t) const override;

    double phi(const double& t) const override;
    array_t phi(const array_t& t) const override;

    double pr(const double& t1, const double& t2, const double& eff_poly) const override;
    array_t pr(const array_t& t1, const array_t& t2, const array_t& eff_poly) const override;

    double eff_poly(const double& p1, const double& t1, const double& p2, const double& t2) const override;
    array_t eff_poly(const array_t& p1, const array_t& t1, const array_t& p2, const array_t& t2) const override;

    double r() const override;

    double t_f_pr(const double& pr_,
                  const double& t1_,
                  const double& eff_poly_,
                  double tol,
                  std::size_t max_iter = 30) const override;

    double t_f_h(const double& h_, double tol, std::size_t max_iter = 30) const override;

    double t_f_phi(const double& phi_, double tol, std::size_t max_iter = 30) const override;

    double static_t(const double& tt_, const double& mach_, double tol, std::size_t max_iter = 30) const override;

    double mach_f_wqa(
        const double& pt, const double& tt, const double& wqa, double tol, std::size_t max_iter = 30) const override;
};

double
PyThermo::gamma(const double& t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface<double>, gamma, t);
}

array_t
PyThermo::gamma(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface<array_t>, v_gamma, t);
}

double
PyThermo::cp(const double& t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface<double>, cp, t);
}

array_t
PyThermo::cp(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface<array_t>, v_cp, t);
}

double
PyThermo::h(const double& t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface<double>, h, t);
}

array_t
PyThermo::h(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface<array_t>, v_h, t);
}

double
PyThermo::phi(const double& t) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface<double>, phi, t);
}

array_t
PyThermo::phi(const array_t& t) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface<array_t>, v_phi, t);
}

double
PyThermo::pr(const double& t1, const double& t2, const double& eff_poly) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface<double>, pr, t1, t2, eff_poly);
}

array_t
PyThermo::pr(const array_t& t1, const array_t& t2, const array_t& eff_poly) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface<array_t>, v_pr, t1, t2, eff_poly);
}

double
PyThermo::eff_poly(const double& p1, const double& t1, const double& p2, const double& t2) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface<double>, eff_poly, p1, t1, p2, t2);
}

array_t
PyThermo::eff_poly(const array_t& p1, const array_t& t1, const array_t& p2, const array_t& t2) const
{
    PYBIND11_OVERLOAD_PURE(array_t, ThermoInterface<array_t>, v_eff_poly, p1, t1, p2, t2);
}

double
PyThermo::r() const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoInterface<double>, r, );
}

double
PyThermo::t_f_pr(const double& pr_, const double& t1_, const double& eff_poly_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoExtendedInterface<double>, t_from_pr, pr_, t1_, eff_poly_, tol, max_iter);
}

double
PyThermo::t_f_h(const double& h_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoExtendedInterface<double>, t_from_h, h_, tol, max_iter);
}

double
PyThermo::t_f_phi(const double& phi_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoExtendedInterface<double>, t_from_phi, phi_, tol, max_iter);
}

double
PyThermo::static_t(const double& tt_, const double& mach_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoExtendedInterface<double>, static_t, tt_, mach_, tol, max_iter);
}

double
PyThermo::mach_f_wqa(const double& pt_, const double& tt_, const double& wqa_, double tol, std::size_t max_iter) const
{
    PYBIND11_OVERLOAD_PURE(double, ThermoExtendedInterface<double>, mach_f_wqa, pt_, tt_, wqa_, tol, max_iter);
}

template class thermo::ThermoExtendedInterface<double>;
template class thermo::ThermoInterface<array_t>;


void
thermo_base(py::module_& m)
{
    using namespace py::literals;

    py::class_<thermo::ThermoExtendedInterface<double>, std::shared_ptr<thermo::ThermoExtendedInterface<double>>>(
        m, "ThermoExtendedInterfaceDouble");


    py::class_<thermo::ThermoInterface<array_t>, std::shared_ptr<thermo::ThermoInterface<array_t>>>(
        m, "ThermoInterfaceArr");

    //  Binding of the abstract/interface class PyThermoInterface
    auto g = py::class_<PyThermoInterface,
                        thermo::ThermoExtendedInterface<double>,
                        thermo::ThermoInterface<array_t>,
                        PyThermo,
                        std::shared_ptr<PyThermoInterface>>(m, "Thermo");
    g.def(py::init<>());

    g.def_property_readonly("constant", &ThermoInterface<double>::r, "Gas constant");
    g.def_property_readonly("r", &ThermoInterface<double>::r, "Gas constant");

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

    py::register_exception<convergence_error>(m, "ConvergenceError", PyExc_RuntimeError);
    py::register_exception<domain_error>(m, "DomainError", PyExc_RuntimeError);
}
