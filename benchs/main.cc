#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"

#include <chrono>
#include <functional>
#include <iostream>


namespace thermo
{
    template <class G>
    void benchmark_scalar(G&& gas, std::size_t size, long ntimes)
    {
        double res;
        srand(time(nullptr));
        double t1 = rand() % 30 + 273.15;
        double t2 = rand() % 30 + 350.;
        double p1 = rand() % 30 + 101325.;
        double p2 = rand() % 30 + 500000.;
        double h = rand() % 30 + 273000.;
        double phi = rand() % (50 + 1) + 6000.;
        double pr = (rand() % 10) / 1000. + 0.75;
        double wqa = rand() % 30 + 180.;

        auto bench_cp = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.cp(t1);
                }
            };
        };

        auto bench_gamma = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.gamma(t1);
                }
            };
        };

        auto bench_h = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.h(t1);
                }
            };
        };

        auto bench_phi = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.phi(t1);
                }
            };
        };

        auto bench_r = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.r();
                }
            };
        };

        auto bench_pr = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.pr(t1, t2, 0.82);
                }
            };
        };

        auto bench_t_from_h = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.t_f_h(h, 1e-8);
                }
            };
        };

        auto bench_t_from_phi = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.t_f_phi(phi, 1e-8);
                }
            };
        };

        auto bench_t_from_pr = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.t_f_pr(pr, t1, 0.82, 1e-8);
                }
            };
        };

        auto bench_eff_poly = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.eff_poly(p1, t1, p2, t2);
                }
            };
        };

        auto bench_pout = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t j = 0; j < size; ++j)
                {
                    res = gas.pr(t1, t2, gas.eff_poly(p1, t1, p2, t2)) * p1;
                }
            };
        };

        auto bench_mach_f_wqa = [&](const long& ntimes) -> void
        {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t j = 0; j < size; ++j)
                {
                    res = gas.mach_f_wqa(p1, t1, wqa, 1e-8, 10);
                }
            };
        };

        auto timeit = [&](std::function<void(long)> f, long ntimes) -> double
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            f(ntimes);
            auto t2 = std::chrono::high_resolution_clock::now();
            return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
                   / static_cast<double>(ntimes * size);
        };

        std::cout << "cp -> " << timeit(bench_cp, ntimes) << " ns" << std::endl;
        std::cout << "gamma -> " << timeit(bench_gamma, ntimes) << " ns" << std::endl;
        std::cout << "h -> " << timeit(bench_h, ntimes) << " ns" << std::endl;
        std::cout << "phi -> " << timeit(bench_phi, ntimes) << " ns" << std::endl;
        std::cout << "r -> " << timeit(bench_r, ntimes) << " ns" << std::endl;
        std::cout << "pr -> " << timeit(bench_pr, ntimes) << " ns" << std::endl;
        std::cout << "t_f_h -> " << timeit(bench_t_from_h, ntimes) << " ns" << std::endl;
        std::cout << "t_f_phi -> " << timeit(bench_t_from_phi, ntimes) << " ns" << std::endl;
        std::cout << "t_f_pr -> " << timeit(bench_t_from_pr, ntimes) << " ns" << std::endl;
        std::cout << "eff_poly -> " << timeit(bench_eff_poly, ntimes) << " ns" << std::endl;
        std::cout << "mach_f_wqa -> " << timeit(bench_mach_f_wqa, ntimes / 2) << " ns" << std::endl;
        std::cout << "p_out -> " << timeit(bench_pout, ntimes) << " ns" << std::endl;
    }
}  // namespace thermo

int
main()
{
    using namespace thermo;
    {
        IdealGas gas(287.05287, 1004.685045);
        std::cout << "\nIdealGas Scalar benchmark" << std::endl;
        benchmark_scalar(gas, 100000, 50);
    }

    {
        PolyGas gas(PolyGasProps<double, 8>({ 0., 0., 0., 0., 0., 0., 0., 1004.4 }, 0., 287.05287));
        std::cout << "\nPolyGas Scalar benchmark" << std::endl;
        benchmark_scalar(gas, 100000, 50);
    }

    return 0;
}
