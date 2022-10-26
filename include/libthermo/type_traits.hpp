// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_TYPE_TRAITS_HPP
#define LIBTHERMO_TYPE_TRAITS_HPP

#include <array>

// clang-format off
#ifdef LIBTHERMO_USE_XTENSOR
    #include "xtensor/xcontainer.hpp"
    #include "xtensor/xtensor.hpp"
    #include "xtensor/xfunction.hpp"
    #include "xtensor/xexpression.hpp"
    #include "xsimd/xsimd.hpp"

    #define IS_NOT_XEXPRESSION std::enable_if_t<!xt::is_xexpression<E>::value, int> = 0
    #define IS_XEXPRESSION std::enable_if_t<xt::is_xexpression<E>::value, int> = 0

    #define IS_NOT_XTENSOR std::enable_if_t<!xt::detail::is_container<T>::value, int> = 0
    #define IS_XTENSOR std::enable_if_t<xt::detail::is_container<T>::value, int> = 0

    #define IS_NOT_XTENSOR_ std::enable_if_t<!xt::detail::is_container<T>::value, int>
    #define IS_XTENSOR_ std::enable_if_t<xt::detail::is_container<T>::value, int>
#else
    #include <type_traits>

    #define IS_NOT_XTENSOR std::enable_if_t<std::true_type(), int> = 0
    #define IS_NOT_XTENSOR_ std::enable_if_t<std::true_type(), int>

    #define IS_XTENSOR std::enable_if_t<std::false_type(), int> = 0
    #define IS_XTENSOR_ std::enable_if_t<std::false_type(), int>

    #define IS_NOT_XEXPRESSION std::enable_if_t<std::true_type(), int> = 0
    #define IS_XEXPRESSION std::enable_if_t<std::false_type(), int> = 0
#endif
// clang-format on

namespace thermo
{
    template <class T, std::size_t N>
    auto array_size_impl(const std::array<T, N>&) -> std::integral_constant<std::size_t, N>;

    template <class Array>
    using array_size = decltype(array_size_impl(std::declval<const Array&>()));

    template <class Array>
    constexpr auto static_size() -> decltype(array_size<Array>::value)
    {
        return array_size<Array>::value;
    }

    template <typename T, typename = void>
    struct is_iterator
    {
        static constexpr bool value = false;
    };

    template <typename T>
    struct is_iterator<
        T,
        typename std::enable_if<
            !std::is_same<typename std::iterator_traits<T>::value_type, void>::value>::type>
    {
        static constexpr bool value = true;
    };
}

#endif
