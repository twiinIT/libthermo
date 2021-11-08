# libthermo

A _fast_ C++ thermodynamical gas library.

`libthermo` targets high-performance computing using really simple/simplistic thermodynamical modelings instead of complex and more accurate approaches (such as https://cantera.org/).

It's well fitted for industrial simulation where the base blocks need to be efficient to not become bottlenecks.

## Modelings

### Properties

The library provides modelings for the following properties:

|  Short  |          Long          |
| :-----: | :--------------------: |
|   Cp    | Specific heat pressure |
|  Gamma  |  Specific heat ratio   |
|    r    |      Gas constant      |
|    H    |        Enthalpy        |
|   Phi   |        Entropy         |
|   PR    |     Pressure ratio     |
| EffPoly | Polytropic efficiency  |

### Properties

For now, only ideal gas is implemented in `IdealGas` class.

It has constant specific heat pressure, thus ratio.

## Benchmarks

Benchmarks can be ran using:

- C++: `thermo_benchs` CMake target
- Python: `ipython python_benchs.ipy`

The metrics are computed over 1M elements for operation only (allocation is not timed), and are given in nanoseconds/element:

### IdealGas

|  Case   | Pure NumPy | C++ loop on vector | C++ xtensor | + xsimd | + tbb |     | pythermo\* |
| :-----: | :--------: | :----------------: | :---------: | :-----: | :---: | --- | :--------: |
|   Cp    |    0.27    |        0.12        |    0.32     |  0.12   | 0.09  |     |            |
|  Gamma  |    0.30    |        0.12        |    0.31     |  0.12   | 0.09  |     |            |
|    r    |    0.28    |        0.13        |             |         | 0.11  |     |            |
|    H    |    0.73    |        0.45        |    0.80     |  0.75   | 0.24  |     |            |
|   Phi   |    4.54    |        4.99        |    4.35     |  2.97   | 0.31  |     |    0.30    |
|   PR    |    11.5    |       14.41        |    13.65    |  8.28   | 0.61  |     |    0.53    |
| EffPoly |    17.1    |       12.18        |    14.87    |  9.87   | 1.02  |     |    2.68    |

1k elements:
|  Case   | Pure NumPy | C++ loop on vector | C++ xtensor | + xsimd | + tbb |     | pythermo\* |
| :-----: | :--------: | :----------------: | :---------: | :-----: | :---: | --- | :--------: |
|   Cp    |            |        0.12        |    0.09     |  0.12   | 0.12  |     |            |
|  Gamma  |            |        0.12        |    0.09     |  0.12   | 0.12  |     |            |
|    r    |            |        0.12        |    0.10     |         |       |     |            |
|    H    |            |        0.13        |    0.09     |  0.14   | 10.2  |     |            |
|   Phi   |            |        5.01        |    4.16     |  3.07   | 11.4  |     |    0.30    |
|   PR    |            |       13.67        |    13.35    |  8.66   | 12.3  |     |    0.53    |
| EffPoly |            |       11.47        |    10.16    |  7.12   | 12.3  |     |    2.68    |


\*pythermo = Python bindings, incl. `xsimd` and `tbb` (vectorization and multithreading)

Conclusions are:

- Pure `NumPy` (vectorized), C++ loop over a 1D buffer allocated on the heap (std::vector) and `xtensor` are directly comparable
  - `xtensor` shows performance issues on very simples cases (`Cp`, `gamma`, `H`): to be investigated
- `xsimd` and `tbb`, turned on by simply adding 2 `CMake` flags give a huge benefit
  - up to x1.65 speed-up for `xsimd`
  - up to x15 speed-up for `tbb`, on a 12-cores/24-threads AMD Ryzen 3900XT
- using `xtensor-python` Python bindings to operate on `NumPy` buffers gives the same performance as pure C++ code

### PolyGas

1M elements:
|  Case   | Pure NumPy | C++ loop on vector | C++ xtensor | + xsimd | + tbb |     | pythermo\* |
| :-----: | :--------: | :----------------: | :---------: | :-----: | :---: | --- | :--------: |
|   Cp    |    0.27    |        0.50        |    0.32     |  0.32   | 0.25  |     |            |
|  Gamma  |    0.30    |        0.60        |    0.31     |  0.31   | 1.96  |     |            |
|    r    |    0.28    |        0.09        |             |         |       |     |            |
|    H    |    0.73    |        0.51        |    0.80     |  0.75   | 0.26  |     |            |
|   Phi   |    4.54    |        7.99        |    4.35     |  2.97   | 0.37  |     |    0.30    |
|   PR    |    11.5    |       21.98        |    13.65    |  8.28   | 0.64  |     |    0.53    |
| EffPoly |    17.1    |       19.49        |    14.87    |  9.87   | 1.03  |     |    2.68    |

1k elements:
|  Case   | Pure NumPy | C++ loop on vector | C++ xtensor | + xsimd | + tbb |     | pythermo\* |
| :-----: | :--------: | :----------------: | :---------: | :-----: | :---: | --- | :--------: |
|   Cp    |    0.27    |        0.29        |    2.31     |  1.26   | 11.2  |     |            |
|  Gamma  |    0.30    |        0.51        |    3.09     |  1.67   | 22.7  |     |            |
|    r    |    0.28    |        0.12        |             |         |       |     |            |
|    H    |    0.73    |        0.36        |    2.82     |  1.56   | 11.5  |     |            |
|   Phi   |    4.54    |        7.84        |    6.15     |  3.62   | 12.3  |     |    0.30    |
|   PR    |    11.5    |       21.97        |    24.34    |  10.5   | 15.2  |     |    0.53    |
| EffPoly |    17.1    |       19.45        |    18.06    |  10.2   | 14.9  |     |    2.68    |
