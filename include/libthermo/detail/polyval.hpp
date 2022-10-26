// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_DETAIL_POLYVAL_HPP
#define LIBTHERMO_DETAIL_POLYVAL_HPP

#include "libthermo/type_traits.hpp"

#include <vector>


namespace
{
    template <int D>
    struct poly_eval
    {
        template <class E, class It, IS_NOT_XEXPRESSION>
        static inline auto eval(E&& x, It coeff)
        {
            auto c = coeff++;
            return std::fma(x, poly_eval<D - 1>::eval(x, coeff), *c);
        };

#ifdef LIBTHERMO_USE_XTENSOR
        template <class E, class It, IS_XEXPRESSION>
        static inline auto eval(E&& x, It coeff)
        {
            auto c = coeff++;
            return xt::fma(x, poly_eval<D - 1>::eval(x, coeff), *c);
        };
#endif
    };

    template <>
    struct poly_eval<0>
    {
        template <class E, class It>
        static inline double eval(E&&, It coeff)
        {
            return *coeff;
        };
    };

    struct poly_impl
    {
        template <class E,
                  class C,
                  int D = thermo::static_size<C>(),
                  std::enable_if_t<std::is_same_v<std::decay_t<C>, std::array<double, D>>, int> = 0>
        static inline auto eval(E&& x, C&& coeffs)
        {
            return poly_eval<D - 1>::eval(x, coeffs.rbegin());
        }

        template <
            class E,
            class C,
            int D = thermo::static_size<C>(),
            std::enable_if_t<!std::is_same_v<std::decay_t<C>, std::array<double, D>>, int> = 0>
        static inline auto eval(E&& x, C&& coeffs)
        {
            return poly_eval<D - 1>::eval(x, coeffs);
        }

        template <class E, class C>
        static inline double eval(E&& x,
                                  C&& coeffs,
                                  std::enable_if_t<std::is_same_v<C, std::vector<double>>, int> = 0)
        {
            std::remove_const_t<std::remove_reference_t<E>> acc = 0;
            for (auto it = coeffs.rbegin(); it < coeffs.rend(); ++it)
                acc = acc * x + (*it);

            return acc;
        }
    };

}

namespace thermo::detail
{
    template <class E, class C>
    static inline auto polyval(E&& x, C&& coeffs)
    {
        return poly_impl::eval(x, coeffs);
    };

    template <int D, class E, class C>
    static inline auto polyval(E&& x, C&& coeffs)
    {
        return poly_impl::eval<E, C, D>(x, coeffs);
    };
}

#endif
