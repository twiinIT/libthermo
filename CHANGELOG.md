# History

## 0.1.2

- improve `mach_f_wqa` iterative function robustness and performance
- fix CMake config to avoid using selection of inappropriate instruction set for SIMD

## 0.1.1

- fix SFINAE error preventing to build on gcc 12 or MSVC

## 0.1.0

- add `Thermo` interface
- add `IdealGas` and `PolyGas` implementations
- support both scalar and vector operations on most of the calculations
- implement efficient `polyval` method for polynomial evaluations
