set(thermo_benchs_srcs main.cc)

add_executable(thermo_benchs ${thermo_benchs_srcs})

target_compile_options(thermo_benchs PUBLIC "-march=native")

target_link_libraries(thermo_benchs
    PUBLIC
        libthermo
)

target_compile_features(thermo_benchs INTERFACE cxx_std_20)

install(TARGETS thermo_benchs DESTINATION ${BIN_DIR})
