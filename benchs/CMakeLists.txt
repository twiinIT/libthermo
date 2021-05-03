set(thermo_benchs_srcs main.cc)

add_executable(thermo_benchs ${thermo_benchs_srcs})
target_link_libraries(thermo_benchs PUBLIC libthermo)
install(TARGETS thermo_benchs DESTINATION ${BIN_DIR})
