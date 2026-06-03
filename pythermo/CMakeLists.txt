# Copyright (c) Adrien DELSALLE
#
# Distributed under the terms of the BSD 3-Clause License.
#
# The full license is in the file LICENSE, distributed with this software.

cmake_minimum_required(VERSION 3.28)

project(pythermo)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

find_package(libthermo REQUIRED)

find_package(Python 3.8
  REQUIRED COMPONENTS Interpreter Development.Module
  OPTIONAL_COMPONENTS Development.SABIModule
)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

set(PYTHERMO_SRCS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

nanobind_add_module(pythermo_core
    STABLE_ABI LTO
    ${PYTHERMO_SRCS_DIR}/main.cpp
    ${PYTHERMO_SRCS_DIR}/thermo.cpp
    ${PYTHERMO_SRCS_DIR}/ideal_gas.cpp
    ${PYTHERMO_SRCS_DIR}/poly_gas.cpp
)

target_compile_features(pythermo_core PUBLIC cxx_std_20)

target_link_libraries(pythermo_core
    PUBLIC
        libthermo
)

install(TARGETS pythermo_core
        LIBRARY DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/pythermo/)
