#include "libthermo/ideal_gas.h"

#include <iostream>
#include <chrono>
#include <functional>
#include <string>
#include <memory>


namespace libthermo
{
    template <class G>
    void Benchmark(G&& gas, long ntimes)
    {
        double res;

        auto bench_Cp = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.Cp(950.); }; };

        auto bench_Gamma = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.Gamma(950.); }; };

        auto bench_H = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.H(950.); }; };

        auto bench_Phi = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.Phi(950.); }; };

        auto bench_R = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.template R<double>(); }; };

        auto bench_PR = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.PR(288.15, 750., 0.82); }; };

        auto bench_EffPoly = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.EffPoly(101325., 288.15, 5e5, 480.); }; };

        auto timeit = [&](std::function<void(long)> f, long ntimes) -> double
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                f(ntimes);
                auto t2 = std::chrono::high_resolution_clock::now();
                return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count()) / static_cast<double>(ntimes);
            };

        auto name = gas.Name();
        std::cout << "Cp -> " << timeit(bench_Cp, ntimes) << " ns" << std::endl;
        std::cout << "Gamma -> " << timeit(bench_Gamma, ntimes) << " ns" << std::endl;
        std::cout << "H -> " << timeit(bench_H, ntimes) << " ns" << std::endl;
        std::cout << "Phi -> " << timeit(bench_Phi, ntimes) << " ns" << std::endl;
        std::cout << "R -> " << timeit(bench_R, ntimes) << " ns" << std::endl;
        std::cout << "PR -> " << timeit(bench_PR, ntimes) << " ns" << std::endl;
        std::cout << "EffPoly -> " << timeit(bench_EffPoly, ntimes) << " ns" << std::endl;
    }
}

int main()
{
    using namespace libthermo;

    IdealGas gas(287.05287, 1004.685045);
    Benchmark(gas, 200000000);
    return 0;
}
