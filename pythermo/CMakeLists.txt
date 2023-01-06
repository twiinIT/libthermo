# Copyright (c) Adrien DELSALLE
#
# Distributed under the terms of the BSD 3-Clause License.
#
# The full license is in the file LICENSE, distributed with this software.

project(pythermo)

set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_MODULE_PATH})

find_package(xtensor-python REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(TBB REQUIRED)
find_package(NumPy REQUIRED)
find_package(Python COMPONENTS Interpreter Development)

set(PYTHERMO_SRCS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

pybind11_add_module(pythermo_core
    ${PYTHERMO_SRCS_DIR}/main.cpp
    ${PYTHERMO_SRCS_DIR}/thermo.cpp
    ${PYTHERMO_SRCS_DIR}/ideal_gas.cpp
    ${PYTHERMO_SRCS_DIR}/poly_gas.cpp
)

target_link_libraries(pythermo_core
    PUBLIC
        libthermo
    PRIVATE
        xtensor-python
        TBB::tbb
)
target_include_directories(pythermo_core PUBLIC ${NUMPY_INCLUDE_DIRS})

install(
    FILES ${PYTHERMO_SRCS_DIR}/thermo.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/pythermo
)

if (WIN32)
    target_compile_definitions(pythermo_core PUBLIC NOMINMAX)
endif ()

install(TARGETS pythermo_core
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pythermo/)
