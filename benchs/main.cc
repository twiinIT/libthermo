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
    nlohmann::json timeit(std::function<void()> f,
                          std::size_t size,
                          unsigned long repeat,
                          unsigned long number)
    {
        xt::xtensor<double, 1>::shape_type times_shape = { number };
        xt::xtensor<double, 1> times(times_shape, 1.);

        for (unsigned long i = 0; i < number; ++i)
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

    nlohmann::json run_benchmark(const std::map<std::string, std::function<void()>>& map,
                                 std::size_t size,
                                 unsigned long repeat,
                                 unsigned long number)
    {
        nlohmann::json j;
        j["N"] = size;
        for (auto& [name, f] : map)
        {
            j[name] = timeit(f, size, repeat, number);
            std::cout << name << " --> " << j[name] << " ns" << std::endl;
        }
        return j;
    }

    void write_results(const std::string& gas_name, const nlohmann::json& j)
    {
        fs::path results_f = "benchs/results/benchs_results.json";
        nlohmann::json existing_j;
        if (fs::exists(results_f))
        {
            std::ifstream current_file(results_f);
            current_file >> existing_j;
        }

        nlohmann::json results;
        results[gas_name] = merge_results(existing_j[gas_name], j);

        fs::create_directories(results_f.parent_path());
        std::ofstream out_file(results_f);
        out_file << results.dump(2);
        out_file.close();
    }

    template <class G>
    void benchmark_single_value(G&& gas, unsigned long repeat, unsigned long number)
    {
        double res;
        double t1 = rand() % 30 + 273.15;
        double t2 = rand() % 30 + 350.;
        double p1 = rand() % 30 + 101325.;
        double p2 = rand() % 30 + 500000.;

        auto bench_Cp = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.Cp(t1);
            };
        };

        auto bench_Gamma = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.Gamma(t1);
            };
        };

        auto bench_H = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.H(t1);
            };
        };

        auto bench_Phi = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.Phi(t1);
            };
        };

        auto bench_R = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.template R<double>();
            };
        };

        auto bench_PR = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.PR(t1, t2, 0.82);
            };
        };

        auto bench_EffPoly = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.EffPoly(p1, t1, p2, t2);
            };
        };

        auto bench_Pout = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                res = gas.PR(t1, t2, gas.EffPoly(p1, t1, p2, t2)) * p1;
            };
        };

        std::map<std::string, std::function<void()>> map
            = { { "Cp", bench_Cp },    { "Gamma", bench_Gamma }, { "H", bench_H },
                { "Phi", bench_Phi },  { "PR", bench_PR },       { "EffPoly", bench_EffPoly },
                { "Pout", bench_Pout } };

        auto j = run_benchmark(map, 1, repeat, number);
        write_results(gas.name(), j);
    }

    template <class G>
    void benchmark_loop(G&& gas, std::size_t size, unsigned long repeat, unsigned long number)
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

        auto bench_Pout = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.PR(t1[i], t2[i], gas.EffPoly(p1[i], t1[i], p2[i], t2[i])) * p1[i];
                }
            };
        };

        std::map<std::string, std::function<void()>> map
            = { { "Cp", bench_Cp },    { "Gamma", bench_Gamma }, { "H", bench_H },
                { "Phi", bench_Phi },  { "PR", bench_PR },       { "EffPoly", bench_EffPoly },
                { "Pout", bench_Pout } };

        auto j = run_benchmark(map, size, repeat, number);

        write_results(gas.name(), j);
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

        auto bench_Pout = [&]() -> void {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.PR(t1, t2, gas.EffPoly(p1, t1, p2, t2)) * p2;
            };
        };

        std::map<std::string, std::function<void()>> map
            = { { "Cp", bench_Cp },    { "Gamma", bench_Gamma }, { "H", bench_H },
                { "Phi", bench_Phi },  { "PR", bench_PR },       { "EffPoly", bench_EffPoly },
                { "Pout", bench_Pout } };

        auto j = run_benchmark(map, size, repeat, number);

        write_results(gas.name(), j);
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

        std::cout << "\nSingle value tests" << std::endl;
        benchmark_single_value(gas, 1000000, 100);

        std::cout << "\nLoop on scalar vector tests" << std::endl;
        benchmark_loop(gas, 1000000, 10, 100);
        // benchmark_loop(gas, 100, 1000000);

#ifdef LIBTHERMO_USE_XTENSOR
        std::cout << "\nVector tests 1M" << std::endl;
        benchmark_vector(gas, 1000000, 10, 100);

        std::cout << "\nVector tests 1k" << std::endl;
        benchmark_vector(gas, 10000, 10, 500);
#endif

#undef TIME_F
#undef TIME
    }
    // std::cout << XSIMD_X86_INSTR_SET << std::endl;

    return 0;
}
