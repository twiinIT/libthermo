# Copyright (c) Adrien DELSALLE
#
# Distributed under the terms of the BSD 3-Clause License.
#
# The full license is in the file LICENSE, distributed with this software.

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})


find_package(xtensor-python REQUIRED)
message(STATUS "Found xtensor-python: ${xtensor-python_INCLUDE_DIRS}/xtensor-python")

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

find_package(NumPy REQUIRED)
message(STATUS "Found NumPy: v${NUMPY_VERSION} ${NUMPY_INCLUDE_DIRS}")


pybind11_add_module(pythermo ideal.cpp)
target_include_directories(pythermo PUBLIC ${NUMPY_INCLUDE_DIRS})
target_link_libraries(pythermo PUBLIC libthermo xtensor-python tbb)

target_compile_options(pythermo PUBLIC "-march=native")
#install(TARGETS pythermo DESTINATION ${SP_DIR})
