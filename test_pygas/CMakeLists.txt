project(test_pygas)

find_package(pybind11 CONFIG REQUIRED)
message(STATUS "Found pybind11 v${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")

# -- Set the Python site-package installation directory
find_package(PythonInterp)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils.sysconfig import get_python_lib; import os;print(get_python_lib())"
                OUTPUT_VARIABLE PY_SP_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)

set(SP_DIR "${PY_SP_DIR}" CACHE PATH "Site-package directory")
message("-- CMake install prefix: " ${CMAKE_INSTALL_PREFIX})
message("  -> site-packages: " ${SP_DIR})

pybind11_add_module(test_pygas main.cpp)
target_include_directories(test_pygas PUBLIC ${NUMPY_INCLUDE_DIRS})
target_link_libraries(test_pygas PUBLIC libthermo)

target_compile_options(test_pygas PUBLIC "-march=native")
install(TARGETS test_pygas DESTINATION ${SP_DIR})