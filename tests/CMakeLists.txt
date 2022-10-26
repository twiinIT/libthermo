# Copyright (c) Adrien DELSALLE
#
# Distributed under the terms of the BSD 3-Clause License.
#
# The full license is in the file LICENSE, distributed with this software.

find_package(GTest)

include_directories(${GTEST_INCLUDE_DIRS})
link_directories(${GTEST_LIBRARIES})

message(STATUS "gtest include dirs: " ${GTEST_INCLUDE_DIRS})
message(STATUS "gtest lib dirs: " ${GTEST_LIBRARIES})

set(libthermo_tests_srcs 
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ideal_gas.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/test_poly_gas.cc
)
add_executable(test_libthermo ${libthermo_tests_srcs} gtest_main.cpp)
target_link_libraries(test_libthermo PUBLIC libthermo GTest::GTest)
install(TARGETS test_libthermo DESTINATION ${BIN_DIR})
