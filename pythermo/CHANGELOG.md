# History

## 0.2.0

- split and template interface between base functions and extended (inverse) ones
- add bindings for all methods of `IdealGas` and `PolyGas8`
- fix headers installation
- fix `IdealGas.mach_f_wqa` iterative function
- fix CMake config files to add Boost dependency

## 0.1.2

- improve `mach_f_wqa` iterative function robustness and performance

## 0.1.1

- fix SFINAE error preventing to build on gcc 12 or MSVC

## 0.1.0

- add Python bindings for `libthermo`, called `pythermo`
- add `IdealGas` bindings
- add `PolyGas<PolyGasProps<8>>` for template specialization of the class providing a 8th order polynomial definition of Cp(s)
