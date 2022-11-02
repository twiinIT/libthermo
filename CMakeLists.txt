# Copyright (c) Adrien DELSALLE
#
# Distributed under the terms of the BSD 3-Clause License.
#
# The full license is in the file LICENSE, distributed with this software.

set(CMAKE_CXX_STANDARD 14)
cmake_minimum_required (VERSION 3.9)
cmake_policy(VERSION 3.9)

project (libthermo)

# Build options
option(BUILD_TESTS "Enable C++ tests for libthermo" OFF)
option(BUILD_BENCHS "Enable C++ benchmarks for libthermo" OFF)
option(BUILD_PY "Build Python bindings of libthermo" OFF)
option(LIBTHERMO_USE_XTENSOR "Enable C++ xtensor for libthermo" OFF)

# Set installation directories (CMAKE_INSTALL_INCLUDEDIR, CMAKE_INSTALL_BINDIR, CMAKE_INSTALL_LIBDIR)
include(GNUInstallDirs)
message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})

set(CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for libthermoConfig.cmake")

message("  -> binaries: " ${CMAKE_INSTALL_BINDIR})
message("  -> libs: " ${CMAKE_INSTALL_LIBDIR})
message("  -> includes: " ${CMAKE_INSTALL_INCLUDEDIR})
message("  -> cmake config: " ${CMAKECONFIG_INSTALL_DIR})

# Add helpers for cmake configuration
include(CMakePackageConfigHelpers)

set(LIBTHERMO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(libthermo_headers
    ${LIBTHERMO_INCLUDE_DIR}/libthermo/thermo.hpp
    ${LIBTHERMO_INCLUDE_DIR}/libthermo/exceptions.hpp
    ${LIBTHERMO_INCLUDE_DIR}/libthermo/ideal_gas.hpp
    ${LIBTHERMO_INCLUDE_DIR}/libthermo/poly_gas.hpp
    ${LIBTHERMO_INCLUDE_DIR}/libthermo/math_utils.hpp
    ${LIBTHERMO_INCLUDE_DIR}/libthermo/version.hpp
)

add_library(libthermo INTERFACE)
target_compile_features(libthermo INTERFACE cxx_std_17)

target_include_directories(
    libthermo
    INTERFACE
        $<BUILD_INTERFACE:${LIBTHERMO_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>
)

find_package(Boost COMPONENTS math)

if (LIBTHERMO_USE_XTENSOR)
    find_package(xtensor REQUIRED)
    target_compile_definitions(libthermo
        INTERFACE
            LIBTHERMO_USE_XTENSOR
    )
    target_link_libraries(libthermo INTERFACE xtensor)
endif ()

target_link_libraries(libthermo
    PUBLIC           
        ${Boost_LIBRARIES}
)

install(TARGETS libthermo 
        EXPORT ${PROJECT_NAME}-targets)

install(FILES ${libthermo_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/libthermo)

if(BUILD_TESTS)
    add_subdirectory(tests)
endif()

if(BUILD_BENCHS)
    add_subdirectory(benchs)
endif()

if(BUILD_PY)
    add_subdirectory(pythermo)
endif()

#add_subdirectory(test_pygas)

# Configuration
# =============
set(ConfigFileInstallDir ${CMAKECONFIG_INSTALL_DIR})

# Definition
configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${ConfigFileInstallDir})

# Installation
install(EXPORT ${PROJECT_NAME}-targets 
        FILE ${PROJECT_NAME}Targets.cmake
        DESTINATION ${ConfigFileInstallDir})

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        DESTINATION ${ConfigFileInstallDir})
