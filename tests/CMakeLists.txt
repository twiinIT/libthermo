# Google Tests #
if(WIN32)
    find_package(GTest)
else()
    add_library(GTest::GTest INTERFACE IMPORTED)
    target_link_libraries(GTest::GTest INTERFACE gtest)
    add_library(GTest::Main INTERFACE IMPORTED)
    target_link_libraries(GTest::Main INTERFACE gtest_main)
endif()

include_directories(${GTEST_INCLUDE_DIRS})
link_directories(${GTEST_LIBRARIES})

message(STATUS "gtest include dirs: " ${GTEST_INCLUDE_DIRS})
message(STATUS "gtest lib dirs: " ${GTEST_LIBRARIES})

set(libthermo_tests_srcs 
    ${CMAKE_CURRENT_SOURCE_DIR}/test_ideal_gas.cc
)
add_executable(test_libthermo ${libthermo_tests_srcs} gtest_main.cpp)
target_link_libraries(test_libthermo PUBLIC libthermo GTest::GTest)
install(TARGETS test_libthermo DESTINATION ${BIN_DIR})
