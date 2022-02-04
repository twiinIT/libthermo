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

# -- Set the Python site-package installation directory
find_package(PythonInterp)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; import os;print(get_python_lib())"
                OUTPUT_VARIABLE PY_SP_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

set(SP_DIR "${PY_SP_DIR}" CACHE PATH "Site-package directory")
message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> site-packages: " ${SP_DIR})

pybind11_add_module(pythermo main.cpp thermo.hpp thermo.cpp ideal_gas.cpp real_gas.cpp)
target_include_directories(pythermo PUBLIC ${NUMPY_INCLUDE_DIRS})
target_link_libraries(pythermo PRIVATE xtensor-python tbb)
target_link_libraries(pythermo PUBLIC libthermo)

target_compile_options(pythermo PUBLIC "-march=native")
install(TARGETS pythermo DESTINATION ${SP_DIR})
