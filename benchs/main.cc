#include "libthermo/ideal_gas.h"
#include "libthermo/real_gas.h"

#ifdef LIBTHERMO_USE_XTENSOR
#include "xtensor/xtensor.hpp"
#endif

#include <iostream>
#include <chrono>
#include <functional>
#include <string>
#include <memory>
#include <vector>

namespace libthermo
{
    template <class G>
    void benchmark_single_value(G&& gas, std::size_t size, long ntimes)
    {
        double res;
        double t1 = rand() % 30 + 273.15;
        double t2 = rand() % 30 + 350.;
        double p1 = rand() % 30 + 101325.;
        double p2 = rand() % 30 + 500000.;

        auto bench_Cp = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res = gas.Cp(t1);
                    }
                };
            };

        auto bench_Gamma = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res = gas.Gamma(t1);
                    }
                };
            };

        auto bench_H = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res = gas.H(t1);
                    }
                };
            };

        auto bench_Phi = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res = gas.Phi(t1);
                    }
                };
            };

        auto bench_R = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res = gas.template R<double>();
                    }
                };
            };

        auto bench_PR = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res = gas.PR(t1, t2, 0.82);
                    }
                };
            };

        auto bench_EffPoly = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res = gas.EffPoly(p1, t1, p2, t2);
                    }
                };
            };

        auto bench_pout = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t j = 0; j < size; ++j)
                    {
                        res = gas.PR(t1, t2, gas.EffPoly(p1, t1, p2, t2)) * p1;
                    }
                };
            };

        auto timeit = [&](std::function<void(long)> f, long ntimes) -> double
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                f(ntimes);
                auto t2 = std::chrono::high_resolution_clock::now();
                return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count()) / static_cast<double>(ntimes * size);
            };

        std::cout << "Cp -> " << timeit(bench_Cp, ntimes*100) << " ns" << std::endl;
        std::cout << "Gamma -> " << timeit(bench_Gamma, ntimes*100) << " ns" << std::endl;
        std::cout << "H -> " << timeit(bench_H, ntimes*10) << " ns" << std::endl;
        std::cout << "Phi -> " << timeit(bench_Phi, ntimes) << " ns" << std::endl;
        std::cout << "R -> " << timeit(bench_R, ntimes*10) << " ns" << std::endl;
        std::cout << "PR -> " << timeit(bench_PR, ntimes) << " ns" << std::endl;
        std::cout << "EffPoly -> " << timeit(bench_EffPoly, ntimes) << " ns" << std::endl;
        std::cout << "Pout -> " << timeit(bench_pout, ntimes) << " ns" << std::endl;
    }

    template <class G>
    void benchmark_loop(G&& gas, std::size_t size, long ntimes)
    {
        std::vector<double> res(size);
        std::vector<double> t1(size, 273.15);
        std::vector<double> t2(size, 350.);
        std::vector<double> p1(size, 101325.);
        std::vector<double> p2(size, 500000.);

        for(auto &t : t1)
            t = rand() % 30 + 273.15;

        auto bench_Cp = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.Cp(t1[i]);
                    }
                };
            };

        auto bench_Gamma = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.Gamma(t1[i]);
                    }
                };
            };

        auto bench_H = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.H(t1[i]);
                    }
                };
            };

        auto bench_Phi = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.Phi(t1[i]);
                    }
                };
            };

        auto bench_R = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.template R<double>();
                    }
                };
            };

        auto bench_PR = [&](const long &ntimes) -> void
            { 
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.PR(t1[i], t2[i], 0.82);
                    }
                };
            };

        auto bench_EffPoly = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.EffPoly(p1[i], t1[i], p2[i], t2[i]);
                    }
                };
            };

        auto bench_pout = [&](const long &ntimes) -> void
            {
                for(long i=0; i<ntimes; i++)
                {
                    for(std::size_t i = 0; i < size; ++i)
                    {
                        res[i] = gas.PR(t1[i], t2[i], gas.EffPoly(p1[i], t1[i], p2[i], t2[i])) * p1[i];
                    }
                };
            };

        auto timeit = [&](std::function<void(long)> f, long ntimes) -> double
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                f(ntimes);
                auto t2 = std::chrono::high_resolution_clock::now();
                return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count()) / static_cast<double>(ntimes * size);
            };

        auto ntimes2 = ntimes / 5;
        std::cout << "Cp -> " << timeit(bench_Cp, ntimes) << " ns" << std::endl;
        std::cout << "Gamma -> " << timeit(bench_Gamma, ntimes) << " ns" << std::endl;
        std::cout << "H -> " << timeit(bench_H, ntimes) << " ns" << std::endl;
        std::cout << "Phi -> " << timeit(bench_Phi, ntimes) << " ns" << std::endl;
        //std::cout << "R -> " << timeit(bench_R, ntimes) << " ns" << std::endl;
        std::cout << "PR -> " << timeit(bench_PR, ntimes2) << " ns" << std::endl;
        std::cout << "EffPoly -> " << timeit(bench_EffPoly, ntimes2) << " ns" << std::endl;
        std::cout << "Pout -> " << timeit(bench_pout, ntimes2) << " ns" << std::endl;
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class G>
    void benchmark_vector(G&& gas, std::size_t size, long ntimes)
    {
        xt::xtensor<double, 1>::shape_type shape = {size};
        xt::xtensor<double, 1> res(shape, 1.);

        xt::xtensor<double, 1> t1(shape, 273.15);
        xt::xtensor<double, 1> t2(shape, 350.);

        xt::xtensor<double, 1> p1(shape, 101325.);
        xt::xtensor<double, 1> p2(shape, 500000.);

        auto bench_Cp = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.Cp(t1); }; };

        auto bench_Gamma = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.Gamma(t1); }; };

        auto bench_H = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.H(t1); }; };

        auto bench_Phi = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.Phi(t1); }; };
/*
        auto bench_R = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.template R<xt::xtensor<double, 1> >(); }; };
*/
        auto bench_PR = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.PR(t1, t2, 0.82); }; };

        auto bench_EffPoly = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.EffPoly(p1, t1, p2, t2); }; };

        auto bench_pout = [&](const long &ntimes) -> void
            { for(long i=0; i<ntimes; i++) { res = gas.PR(t1, t2, gas.EffPoly(p1, t1, p2, t2)) * p2; }; };

        auto timeit = [&](std::function<void(long)> f, long ntimes) -> double
            {
                auto t1 = std::chrono::high_resolution_clock::now();
                f(ntimes);
                auto t2 = std::chrono::high_resolution_clock::now();
                return static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>( t2 - t1 ).count()) / static_cast<double>(ntimes * size);
            };

        auto ntimes2 = ntimes / 5;
        std::cout << "Cp -> " << timeit(bench_Cp, ntimes) << " ns" << std::endl;
        std::cout << "Gamma -> " << timeit(bench_Gamma, ntimes) << " ns" << std::endl;
        std::cout << "H -> " << timeit(bench_H, ntimes) << " ns" << std::endl;
        std::cout << "Phi -> " << timeit(bench_Phi, ntimes) << " ns" << std::endl;
        // std::cout << "R -> " << timeit(bench_R, ntimes) << " ns" << std::endl;
        std::cout << "PR -> " << timeit(bench_PR, ntimes) << " ns" << std::endl;
        std::cout << "EffPoly -> " << timeit(bench_EffPoly, ntimes) << " ns" << std::endl;
        std::cout << "Pout -> " << timeit(bench_pout, ntimes) << " ns" << std::endl;
    }
#endif
}

int main()
{
    using namespace libthermo;
    double Tref = 288.15;
    {
        IdealGas gas(287.05287, 1004.685045);

        std::cout << "\n" << "Reference values" << std::endl;
        std::cout << "Cp -> " << gas.Cp(Tref) << std::endl;
        std::cout << "Gamma -> " << gas.Gamma(Tref) << std::endl;
        std::cout << "H -> " << gas.H(Tref) << std::endl;
        std::cout << "Phi -> " << gas.Phi(Tref) << std::endl;

        std::cout << "\nSingle value tests" << std::endl;
        benchmark_single_value(gas, 1000000, 50);

        std::cout << "\nLoop tests" << std::endl;
        benchmark_loop(gas, 1000000, 500);

        std::cout << "\nVector tests" << std::endl;
        benchmark_vector(gas, 1000000, 500);
    }

    {
        RealGas gas(287.05287);

        std::cout << "\n" << "Reference values" << std::endl;
        std::cout << "Cp -> " << gas.Cp(Tref) << std::endl;
        std::cout << "Gamma -> " << gas.Gamma(Tref) << std::endl;
        std::cout << "H -> " << gas.H(Tref) << std::endl;
        std::cout << "Phi -> " << gas.Phi(Tref) << std::endl;

        std::cout << "\nSingle value tests" << std::endl;
        benchmark_single_value(gas, 1000000, 50);

        std::cout << "\nLoop tests" << std::endl;
        benchmark_loop(gas, 1000000, 500);

        std::cout << "\nVector tests" << std::endl;
        benchmark_vector(gas, 1000000, 500);

    }
    // std::cout << XSIMD_X86_INSTR_SET << std::endl;

    return 0;
}
