// Copyright (c) 2021-2023, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include "thermo.hpp"
#include "libthermo/ideal_gas.hpp"

#include <nanobind/ndarray.h>
// #include <nanobind/stl.h>

#include <numeric>
#include <vector>


namespace nb = nanobind;

using namespace thermo;

template <typename T, typename... Ts>
bool
have_same_dimensions(const T& first, const Ts&... rest)
{
    return ((first.ndim() == rest.ndim()) && ...);
}

template <typename T, typename... Ts>
bool
have_same_shapes(const T& first, const Ts&... rest)
{
    bool are_same = have_same_dimensions(first, rest...);
    if (!are_same)
        return false;

    for (size_t i = 0; i < first.ndim(); ++i)
        are_same &= ((first.shape(i) == rest.shape(i)) && ...);

    return are_same;
}

// nb::ndarray<double>
// as_array(nb::handle h)
// {
//     if (nb::isinstance<nb::ndarray<>>(h))
//     {
//         return nb::cast<nb::ndarray<double>>(h);
//     }
//     else if (nb::isinstance<nb::float_>(h) || nb::isinstance<nb::int_>(h))
//     {
//         double* res = (double*) std::malloc(array.size() * sizeof(double));
//         if (!res)
//             throw std::bad_alloc();


//         //   nb::capsule owner(res, [](void* ptr) noexcept { std::free(ptr); });

//         //   return nb::ndarray<nb::numpy, double>(res, { array.size() }, owner);

//         // Wrap scalar into 0-d ndarray
//         nb::ndarray<double> arr{};  // 0-dim
//         arr.data()[0] = ;
//         return arr;
//     }
//     else
//     {
//         throw std::invalid_argument("Expected float, int or ndarray");
//     }
// }

nb::ndarray<double>
broadcast_input(nb::handle input, std::size_t size, std::size_t ndim, const std::size_t* shapes)
{
    if (nb::isinstance<nb::ndarray<double>>(input))
        return nb::cast<nb::ndarray<double>>(input);

    double val = nb::cast<double>(input);

    double* res = (double*) std::malloc(size * sizeof(double));
    if (!res)
        throw std::bad_alloc();

    std::fill(res, res + size, val);

    nb::capsule owner(res, [](void* ptr) noexcept { std::free(ptr); });

    return nb::ndarray<double>{ res, ndim, shapes, owner };
}


namespace pythermo
{
    void ideal_gas(nb::module_& m)
    {
        using namespace nb::literals;

        auto ideal_gas = nb::class_<IdealGas<double>>(m, "IdealGas");
        ideal_gas.def(nb::init<double, double>());

        register_base_interface<nb::class_<IdealGas<double>>, double>(ideal_gas);


        ideal_gas.def("h",
                      [](IdealGas<double>& gas, const nb::ndarray<double>& array, nb::ndarray<double>& res)
                      {
                          //   double* res = (double*) std::malloc(array.size() * sizeof(double));
                          //   if (!res)
                          //       throw std::bad_alloc();

                          gas.h(res.data(), array.data(), array.size());

                          //   nb::capsule owner(res, [](void* ptr) noexcept { std::free(ptr); });

                          //   return nb::ndarray<nb::numpy, double>(res, { array.size() }, owner);
                      });

        ideal_gas.def(
            "eff_poly",
            [](IdealGas<double>& gas,
               const nb::ndarray<double, nb::device::cpu>& p1,
               const nb::ndarray<double, nb::device::cpu>& t1,
               const nb::ndarray<double, nb::device::cpu>& p2,
               const nb::ndarray<double, nb::device::cpu>& t2,
               nb::ndarray<double>& res)
            {
                if (!have_same_shapes(p1, p2, t1, t2, res))
                    throw std::runtime_error("Incompatible shapes");

                gas.eff_poly(res.data(), p1.data(), t1.data(), p2.data(), t2.data(), p1.size());
            },
            "Polytropic efficiency (vectorized)",
            "p1"_a,
            "t1"_a,
            "p2"_a,
            "t2"_a,
            "result"_a);

        ideal_gas.def(
            "eff_poly",
            [](IdealGas<double>& gas,
               nb::handle p1,
               nb::handle t1,
               nb::handle p2,
               nb::handle t2,
               nb::ndarray<double, nb::device::cpu> res)
            {
                std::vector<nb::handle> handles{ p1, t1, p2, t2 };
                std::vector<nb::ndarray<double>> broadcasted(5);
                for (auto i = 0; i < handles.size(); ++i)
                {
                    try
                    {
                        broadcasted[i]
                            = broadcast_input(handles[i], res.size(), res.ndim(), (const std::size_t*) res.shape_ptr());
                    }
                    catch (...)
                    {
                        throw nb::type_error("Expected float, int or ndarray");
                    }
                }

                gas.eff_poly(res.data(),
                             broadcasted[0].data(),
                             broadcasted[1].data(),
                             broadcasted[2].data(),
                             broadcasted[3].data(),
                             res.size());
            },
            "Polytropic efficiency (vectorized)",
            "p1"_a,
            "t1"_a,
            "p2"_a,
            "t2"_a,
            "result"_a);
        // ideal_gas.def("constant", &IdealGas<double>::r, "Gas constant").def("r", &IdealGas<double>::r, "Gas
        // constant");

        // ideal_gas.def(nb::pickle(
        //     [](const PyIdealGas& g) -> nb::tuple {  // __getstate__
        //         return nb::make_tuple(g.r(), g.cp(288.15));
        //     },
        //     [](nb::tuple t) -> PyIdealGas {  // __setstate__
        //         PyIdealGas g(t[0].cast<double>(), t[1].cast<double>());
        //         return g;
        //     }));

        // Define `IdealGas` overloads of inverse methods not requiring iterations
        // Those overloads are NOT part of the `ThermoExtendedInterface` but are provided for convenience
        // note: base class overloads need to be redefined here to not be masked by derived class ones
        // ideal_gas.def(
        //     "static_t", &IdealGas<double>::static_t, "Static temperature from total and Mach number", "tt"_a,
        //     "mach"_a);
        // ideal_gas.def("static_t",
        //               &ThermoExtendedInterface<double>::static_t,
        //               "Static temperature from total and Mach number",
        //               "tt"_a,
        //               "mach"_a,
        //               "tol"_a,
        //               "max_iter"_a = 30);

        // ideal_gas.def("t_f_h", &PyIdealGas::t_f_h, "Temperature from enthalpy", "h"_a);
        // ideal_gas.def("t_f_h",
        //               &ThermoExtendedInterface<double>::t_f_h,
        //               "Temperature from enthalpy",
        //               "h"_a,
        //               "tol"_a,
        //               "max_iter"_a = 30);

        // ideal_gas.def("t_f_phi", &PyIdealGas::t_f_phi, "Temperature from phi function", "phi"_a);
        // ideal_gas.def("t_f_phi",
        //               &ThermoExtendedInterface<double>::t_f_phi,
        //               "Temperature from phi function",
        //               "h"_a,
        //               "tol"_a,
        //               "max_iter"_a = 30);

        // ideal_gas.def("t_f_pr", &PyIdealGas::t_f_pr, "Temperature from pressure ratio", "pr"_a, "t1"_a,
        // "eff_poly"_a); ideal_gas.def("t_f_pr",
        //               &ThermoExtendedInterface<double>::t_f_pr,
        //               "Temperature from pressure ratio",
        //               "pr"_a,
        //               "t1"_a,
        //               "eff_poly"_a,
        //               "tol"_a,
        //               "max_iter"_a = 30);

        // ideal_gas.def("t_f_pr", &PyIdealGas::t_f_pr, "Temperature from pressure ratio", "pr"_a, "t1"_a,
        // "eff_poly"_a); ideal_gas.def("t_f_pr",
        //               &ThermoExtendedInterface<double>::t_f_pr,
        //               "Temperature from pressure ratio",
        //               "pr"_a,
        //               "t1"_a,
        //               "eff_poly"_a,
        //               "tol"_a,
        //               "max_iter"_a = 30);
    }
}