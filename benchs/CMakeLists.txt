set(thermo_benchs_srcs main.cc)

find_package(nlohmann_json 3.11.2 REQUIRED)
find_package(TBB REQUIRED)

add_executable(thermo_benchs ${thermo_benchs_srcs})

target_compile_options(thermo_benchs PUBLIC "-march=native")

add_definitions(-DXTENSOR_TBB_THRESHOLD=5000)
add_definitions(-DXTENSOR_USE_TBB)
add_definitions(-DXTENSOR_USE_XSIMD)

target_link_libraries(thermo_benchs
    PUBLIC
        libthermo
        TBB::tbb
    PRIVATE
        nlohmann_json::nlohmann_json
)

target_compile_features(thermo_benchs INTERFACE cxx_std_20)

install(TARGETS thermo_benchs DESTINATION ${BIN_DIR})
