// Copyright (c) 2021-2022, twiinIT
//
// Distributed under the terms of the BSD 3-Clause License.
//
// The full license is in the file LICENSE, distributed with this software.

#ifndef LIBTHERMO_VERSION_HPP
#define LIBTHERMO_VERSION_HPP

#define LIBTHERMO_VERSION_MAJOR 0
#define LIBTHERMO_VERSION_MINOR 1
#define LIBTHERMO_VERSION_PATCH 1

#define __LIBTHERMO_STRINGIZE_IMPL(s) #s
#define __LIBTHERMO_STRINGIZE(s) __LIBTHERMO_STRINGIZE_IMPL(s)

#define LIBTHERMO_VERSION                                                                          \
    (LIBTHERMO_VERSION_MAJOR * 10000 + LIBTHERMO_VERSION_MINOR * 100 + LIBTHERMO_VERSION_PATCH)
#define LIBTHERMO_VERSION_STRING                                                                   \
    __LIBTHERMO_STRINGIZE(LIBTHERMO_VERSION_MAJOR)                                                 \
    "." __LIBTHERMO_STRINGIZE(LIBTHERMO_VERSION_MINOR) "." __LIBTHERMO_STRINGIZE(                  \
        LIBTHERMO_VERSION_PATCH)

#endif
