#include "libthermo/ideal_gas.h"
#include "libthermo/real_gas.h"

#ifdef LIBTHERMO_USE_XTENSOR
#include "xtensor/xtensor.hpp"
#include "xtensor/xnoalias.hpp"
#endif

#include "nlohmann/json.hpp"

#include <chrono>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <filesystem>


namespace fs = std::filesystem;

namespace libthermo
{
    nlohmann::json timeit(std::function<void()> f, std::size_t size, long repeat, long number)
    {
        xt::xtensor<double, 1>::shape_type times_shape = { number };
        xt::xtensor<double, 1> times(times_shape, 1.);

        for (long i = 0; i < number; ++i)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            f();
            auto t2 = std::chrono::high_resolution_clock::now();
            auto t = static_cast<double>(
                         std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
                     / static_cast<double>(repeat * size);

            times[i] = t;
        }
        return { { xt::mean(times)(), xt::stddev(times)() } };
    }

    nlohmann::json merge_results(const nlohmann::json& a, const nlohmann::json& b)
    {
        if (a.is_null())
            return b;
        if (b.is_null())
            return a;

        nlohmann::json res = a;
        for (auto& [key, value] : b.items())
        {
            if (a.contains(key))
            {
                if (a[key].is_object() && value.is_object())
                    res[key] = merge_results(a[key], value);
                else if (a[key].is_array())
                    res[key].insert_iterator(res[key].end(), value.begin(), value.end());
                else if (a[key].is_number())
                    res[key] = { a[key], value };
            }
            else
            {
                res[key].push_back(value);
            }
        }
        return res;
    }

    template <class G>
    void benchmark_single_value(G&& gas, std::size_t size, long ntimes)
    {
        double res;
        double t1 = rand() % 30 + 273.15;
        double t2 = rand() % 30 + 350.;
        double p1 = rand() % 30 + 101325.;
        double p2 = rand() % 30 + 500000.;

        auto bench_Cp = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.Cp(t1);
                }
            };
        };

        auto bench_Gamma = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.Gamma(t1);
                }
            };
        };

        auto bench_H = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.H(t1);
                }
            };
        };

        auto bench_Phi = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.Phi(t1);
                }
            };
        };

        auto bench_R = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.template R<double>();
                }
            };
        };

        auto bench_PR = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.PR(t1, t2, 0.82);
                }
            };
        };

        auto bench_EffPoly = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res = gas.EffPoly(p1, t1, p2, t2);
                }
            };
        };

        auto bench_pout = [&](const long& ntimes) -> void {
            for (long i = 0; i < ntimes; i++)
            {
                for (std::size_t j = 0; j < size; ++j)
                {
                    res = gas.PR(t1, t2, gas.EffPoly(p1, t1, p2, t2)) * p1;
                }
            };
        };

        auto timeit = [&](std::function<void(long)> f, long ntimes) -> double {
            auto t1 = std::chrono::high_resolution_clock::now();
            f(ntimes);
            auto t2 = std::chrono::high_resolution_clock::now();
            return static_cast<double>(
                       std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
                   / static_cast<double>(ntimes * size);
        };

        std::cout << "Cp -> " << timeit(bench_Cp, ntimes * 100) << " ns" << std::endl;
        std::cout << "Gamma -> " << timeit(bench_Gamma, ntimes * 100) << " ns" << std::endl;
        std::cout << "H -> " << timeit(bench_H, ntimes * 10) << " ns" << std::endl;
        std::cout << "Phi -> " << timeit(bench_Phi, ntimes) << " ns" << std::endl;
        std::cout << "R -> " << timeit(bench_R, ntimes * 10) << " ns" << std::endl;
        std::cout << "PR -> " << timeit(bench_PR, ntimes) << " ns" << std::endl;
        std::cout << "EffPoly -> " << timeit(bench_EffPoly, ntimes) << " ns" << std::endl;
        std::cout << "Pout -> " << timeit(bench_pout, ntimes) << " ns" << std::endl;
    }

    template <class G>
    void benchmark_loop(G&& gas, std::size_t size, long repeat, long number)
    {
        std::vector<double> res(size);
        std::vector<double> t1(size, 273.15);
        std::vector<double> t2(size, 350.);
        std::vector<double> p1(size, 101325.);
        std::vector<double> p2(size, 500000.);

        for (auto& t : t1)
            t = rand() % 30 + 273.15;

        auto bench_Cp = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.Cp(t1[i]);
                }
            };
        };

        auto bench_Gamma = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.Gamma(t1[i]);
                }
            };
        };

        auto bench_H = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.H(t1[i]);
                }
            };
        };

        auto bench_Phi = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.Phi(t1[i]);
                }
            };
        };

        auto bench_R = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.template R<double>();
                }
            };
        };

        auto bench_PR = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.PR(t1[i], t2[i], 0.82);
                }
            };
        };

        auto bench_EffPoly = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.EffPoly(p1[i], t1[i], p2[i], t2[i]);
                }
            };
        };

        auto bench_pout = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.PR(t1[i], t2[i], gas.EffPoly(p1[i], t1[i], p2[i], t2[i])) * p1[i];
                }
            };
        };

        auto number2 = number / 5;
        nlohmann::json j;
        j["N"] = size;

#define TIME_F(NAME)                                                                               \
    j[#NAME] = timeit(bench_ #NAME, size, repeat, number);                                         \
    std::cout << #NAME << " --> " << j[#NAME] << " ns";


        j["Cp"] = TIME_F(Cp);
        j["Gamma"] = TIME_F(Gamma);
        j["H"] = TIME_F(H);
        j["Phi"] = TIME_F(Phi);
        j["PR"] = TIME_F(PR);
        j["EffPoly"] = TIME_F(EffPoly);
        j["Pout"] = TIME_F(Pout);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class G>
    void benchmark_vector(G&& gas, std::size_t size, long repeat, long number)
    {
        xt::xtensor<double, 1>::shape_type shape = { size };
        xt::xtensor<double, 1> res(shape, 1.);

        xt::xtensor<double, 1> t1(shape, 273.15);
        xt::xtensor<double, 1> t2(shape, 350.);

        xt::xtensor<double, 1> p1(shape, 101325.);
        xt::xtensor<double, 1> p2(shape, 500000.);

        xt::xtensor<double, 1> eff(shape, 0.82);

        auto bench_Cp = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.Cp(t1);
            };
        };

        auto bench_Gamma = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.Gamma(t1);
            };
        };

        auto bench_H = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.H(t1);
            };
        };

        auto bench_Phi = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.Phi(t1);
            };
        };
        /*
            auto bench_dPhi = [&]() -> void {
                for (long i = 0; i < repeat; i++)
                {
                    xt::noalias(res) = gas.dPhi(t1, t2);
                };
            };
        */
        /*
                auto bench_R = [&](const long &ntimes) -> void
                    { for(long i=0; i<ntimes; i++) { res = gas.template
           R<xt::xtensor<double, 1> >(); }; };
        */
        auto bench_PR = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.PR(t1, t2, eff);
            };
        };

        auto bench_EffPoly = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.EffPoly(p1, t1, p2, t2);
            };
        };

        auto bench_pout = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.PR(t1, t2, gas.EffPoly(p1, t1, p2, t2)) * p2;
            };
        };

        nlohmann::json j;
        j["N"] = size;
        j["Cp"] = timeit(bench_Cp, size, repeat, number);
        j["Gamma"] = timeit(bench_Gamma, size, repeat, number);
        j["H"] = timeit(bench_H, size, repeat, number);
        j["Phi"] = timeit(bench_Phi, size, repeat, number);
        j["PR"] = timeit(bench_PR, size, repeat, number);
        j["EffPoly"] = timeit(bench_EffPoly, size, repeat, number);
        j["Pout"] = timeit(bench_pout, size, repeat, number);

        std::cout << j.dump() << std::endl;

        fs::path results_f = "benchs/results/benchs_results.json";
        nlohmann::json existing_j;
        if (fs::exists(results_f))
        {
            std::ifstream current_file(results_f);
            current_file >> existing_j;
        }

        nlohmann::json results;
        results[gas.name()] = merge_results(existing_j[gas.name()], j);

        std::ofstream out_file(results_f);
        out_file << results.dump(2);
        out_file.close();
    }
#endif
}  // namespace libthermo

int
main()
{
    using namespace libthermo;
    double Tref = 288.15;
    {
        IdealGas gas(287.05287, 1004.685045);

        std::cout << "\n"
                  << "Reference values" << std::endl;
        std::cout << "Cp -> " << gas.Cp(Tref) << std::endl;
        std::cout << "Gamma -> " << gas.Gamma(Tref) << std::endl;
        std::cout << "H -> " << gas.H(Tref) << std::endl;
        std::cout << "Phi -> " << gas.Phi(Tref) << std::endl;
        /*
                std::cout << "\nSingle value tests" << std::endl;
                benchmark_single_value(gas, 1000000, 50);

                std::cout << "\nLoop tests" << std::endl;
                benchmark_loop(gas, 1000000, 500);

        #ifdef LIBTHERMO_USE_XTENSOR
                std::cout << "\nVector tests" << std::endl;
                benchmark_vector(gas, 1000000, 5000);
        #endif
        */
    }

    {
        RealGas gas(287.05287);

        std::cout << "\n"
                  << "Reference values" << std::endl;
        std::cout << "Cp -> " << gas.Cp(Tref) << std::endl;
        std::cout << "Gamma -> " << gas.Gamma(Tref) << std::endl;
        std::cout << "H -> " << gas.H(Tref) << std::endl;
        std::cout << "Phi -> " << gas.Phi(Tref) << std::endl;
        /*
                std::cout << "\nSingle value tests" << std::endl;
                benchmark_single_value(gas, 1000000, 50);
        std::cout << "\nLoop tests" << std::endl;
        benchmark_loop(gas, 10, 10000000);
        benchmark_loop(gas, 100, 1000000);
        */

#ifdef LIBTHERMO_USE_XTENSOR
        std::cout << "\nVector tests 1M" << std::endl;
        benchmark_vector(gas, 1000000, 1000, 10);

        std::cout << "\nVector tests 1k" << std::endl;
        benchmark_vector(gas, 10000, 1000, 100);
#endif
    }
    // std::cout << XSIMD_X86_INSTR_SET << std::endl;

    return 0;
}
