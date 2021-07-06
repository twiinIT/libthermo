# libthermo

A *fast* C++ thermodynamical library.

`libthermo` targets high-performance computing using really simple/simplistic thermodynamical modelings instead of complex and more accurate approaches (such as https://cantera.org/).

It's well fitted for industrial simulation where the base blocks need to be efficient to not become bottlenecks.

## Benchmarks

Benchmarks can be ran using:
- C++: `thermo_benchs` CMake target
- Python: `ipython python_benchs.ipy`

The metrics are computed over 1M elements for operation only (allocation is not timed), and are given in nanoseconds/element:

|  Case  | Pure NumPy  | C++ loop on vector | C++ xtensor  | + xsimd  | + tbb |  | pythermo* |
|:------:|:-----------:|:------------------:|:------------:|:--------:|:-----:|--|:---------:|
|   Cp   |     0.27    |      0.12          |    0.32      |   0.32   | 0.31  |  |           |
| Gamma  |     0.30    |      0.12          |    0.31      |   0.31   | 0.31  |  |           |
|   r    |     0.28    |      0.45          |              |          |       |  |           |
|   H    |     0.73    |      0.45          |    0.80      |   0.75   | 0.26  |  |           |
|  Phi   |     4.54    |      4.99          |    4.35      |   2.97   | 0.30  |  |  0.30     |
|   PR   |     11.5    |      14.41         |    13.65     |   8.28   | 0.55  |  |  0.53     |
|EffPoly |     17.1    |      16.18         |    14.87     |   9.87   | 2.59  |  |  2.68     |

*pythermo = Python bindings, incl. `xsimd` and `tbb` (vectorization and multithreading)

Conclusions are:
- Pure `NumPy` (vectorized), C++ loop over a 1D buffer allocated on the heap (std::vector) and `xtensor` are directly comparable
  - `xtensor` shows performance issues on very simples cases (`Cp`, `gamma`, `H`): to be investigated
- `xsimd` and `tbb`, turned on by simply adding 2 `CMake` flags give a huge benefit
  - up to x1.65 speed-up for `xsimd`
  - up to x15 speed-up for `tbb`, on a 12-cores/24-threads AMD Ryzen 3900XT
- using `xtensor-python` Python bindings to operate on `NumPy` buffers gives the same performance as pure C++ code
