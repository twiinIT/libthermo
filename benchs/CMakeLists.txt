set(thermo_benchs_srcs 
    # main.cc
    multiprops_array.cc
)

find_package(fmt REQUIRED)
find_package(nlohmann_json 3.11.2 REQUIRED)
# find_package(Cantera REQUIRED)
# find_package(TBB REQUIRED)

add_executable(thermo_benchs ${thermo_benchs_srcs})

target_compile_options(thermo_benchs PUBLIC /utf-8)

add_definitions(-DXTENSOR_TBB_THRESHOLD=5000)
# add_definitions(-DXTENSOR_USE_TBB)
add_definitions(-DXTENSOR_USE_XSIMD)

target_link_libraries(thermo_benchs
    PUBLIC
        libthermo
        # TBB::tbb
    PRIVATE
        nlohmann_json::nlohmann_json
        cantera_shared
        fmt
)

target_compile_features(thermo_benchs INTERFACE cxx_std_20)

install(TARGETS thermo_benchs DESTINATION ${BIN_DIR})
