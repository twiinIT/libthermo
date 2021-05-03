// Copyright (c) Adrien DELSALLE
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#include<assert.h>
#include<vector>

#include "libthermo/math_utils.h"


double math::Polyval(std::vector<double> coeffs, double value)
{
    assert(coeffs.size() > 0);
    double result = coeffs[0];
    for (int i = 1; i < coeffs.size(); i++)
    {
        result = result * value + coeffs[i];
    }
    return result;
}
