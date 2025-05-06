# Copyright (c) Adrien DELSALLE
#
# Distributed under the terms of the BSD 3-Clause License.
#
# The full license is in the file LICENSE, distributed with this software.

cmake_minimum_required(VERSION 3.15...3.29)

project(pythermo)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT
)

find_package(nanobind CONFIG REQUIRED)
find_package(xsimd REQUIRED)
find_package(libthermo REQUIRED)

set(PYTHERMO_SRCS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

nanobind_add_module(
    pythermo_core

    STABLE_ABI
    NB_STATIC

    ${PYTHERMO_SRCS_DIR}/main.cpp
    ${PYTHERMO_SRCS_DIR}/thermo.cpp
    ${PYTHERMO_SRCS_DIR}/ideal_gas.cpp
    # ${PYTHERMO_SRCS_DIR}/poly_gas.cpp
)

# target_include_directories(
#     pythermo_core
#     PUBLIC $<BUILD_INTERFACE:${CONFIGURABLE_INCLUDE_DIR}> $<INSTALL_INTERFACE:include>
# )


# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

target_link_libraries(pythermo_core
    PUBLIC
        libthermo
    PRIVATE
        xsimd
)

# target_include_directories(pythermo_core PUBLIC ${NUMPY_INCLUDE_DIRS})

if (WIN32)
    target_compile_definitions(pythermo_core PUBLIC NOMINMAX)
endif ()

install(TARGETS pythermo_core
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pythermo/)
