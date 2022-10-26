// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_MATH_UTILS_HPP
#define LIBTHERMO_MATH_UTILS_HPP

#include <cmath>
#include <iostream>
#include <vector>


namespace math
{
    constexpr double pi()
    {
        return 3.141592653589793;
    }

    inline double radians(double angle)
    {
        return pi() * angle / 180.;
    }
    inline double degrees(double angle)
    {
        return angle * 180. / pi();
    }

    inline double rpm(double omega)
    {
        return omega * 30. / pi();
    }
    inline double radsec(double omega)
    {
        return omega * pi() / 30.;
    }

    inline double sign(double a)
    {
        if (a >= 0.)
        {
            return 1.;
        }
        else
        {
            return -1.;
        }
    }
    inline double square(double a)
    {
        return a * a;
    }
}
#endif
