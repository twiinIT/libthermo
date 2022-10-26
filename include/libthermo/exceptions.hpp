// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_EXCEPTIONS_HPP
#define LIBTHERMO_EXCEPTIONS_HPP

#include <exception>
#include <string>


namespace thermo
{
    /**
     * Convergence exception
     *
     * This exception is thrown when a numerical algorithm
     * didn't converged at the precision required by the caller.
     */
    class convergence_error : public std::exception
    {
    public:
        inline const char* what() const noexcept
        {
            return "Algorithm did not converged at the specified tolerance";
        };
    };

}

#endif
