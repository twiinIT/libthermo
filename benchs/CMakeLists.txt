set(thermo_benchs_srcs main.cc)

find_package(nlohmann_json 3.11.2 REQUIRED)

add_executable(thermo_benchs ${thermo_benchs_srcs})

target_compile_options(thermo_benchs PUBLIC "-march=native")

add_definitions(-DXTENSOR_TBB_THRESHOLD=5000)

target_link_libraries(thermo_benchs
    PUBLIC
        libthermo
        TBB::tbb
    PRIVATE
        nlohmann_json::nlohmann_json
)
install(TARGETS thermo_benchs DESTINATION ${BIN_DIR})
