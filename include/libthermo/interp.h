
#include "libthermo/polynomial_gas.h"

#include <array>


struct LinearFunc
{
    const std::array<double, 2> coeffs;

    template<class T>
    auto operator()(T&& x) const;
};


template<class T>
auto LinearFunc::operator()(T&& x) const
{
    return std::move(libthermo::polyval<1>(x, coeffs.rbegin()));
}
