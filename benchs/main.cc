#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"

#ifdef LIBTHERMO_USE_XTENSOR  // clang-format off
    #include "xtensor/xtensor.hpp"
    #include "xtensor/xnoalias.hpp"
#endif  // clang-format on

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

namespace thermo
{
    nlohmann::json timeit(std::function<void()> f, std::size_t size, unsigned long repeat, unsigned long number)
    {
        xt::xtensor<double, 1>::shape_type times_shape = { number };
        xt::xtensor<double, 1> times(times_shape, 1.);

        for (long i = 0; i < number; ++i)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            f();
            auto t2 = std::chrono::high_resolution_clock::now();
            auto t = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
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
                    res = gas.mach_f_wqa(p1, t1, wqa, 1e-8);
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

        auto bench_cp = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.cp(t1[i]);
                }
            };
        };

        auto bench_gamma = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.gamma(t1[i]);
                }
            };
        };

        auto bench_h = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.h(t1[i]);
                }
            };
        };

        auto bench_phi = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.phi(t1[i]);
                }
            };
        };

        auto bench_r = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.r();
                }
            };
        };

        auto bench_pr = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.pr(t1[i], t2[i], 0.82);
                }
            };
        };

        auto bench_eff_poly = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.eff_poly(p1[i], t1[i], p2[i], t2[i]);
                }
            };
        };

        auto bench_Pout = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                for (std::size_t i = 0; i < size; ++i)
                {
                    res[i] = gas.pr(t1[i], t2[i], gas.eff_poly(p1[i], t1[i], p2[i], t2[i])) * p1[i];
                }
            };
        };

        auto number2 = number / 5;
        nlohmann::json j;
        j["N"] = size;

#define TIME_F(NAME)                                                                                                   \
    j[#NAME] = timeit(bench_##NAME, size, repeat, number);                                                             \
    std::cout << #NAME << " --> " << j[#NAME] << " ns";


        j["cp"] = TIME_F(cp);
        j["gamma"] = TIME_F(gamma);
        j["h"] = TIME_F(h);
        j["phi"] = TIME_F(phi);
        j["pr"] = TIME_F(pr);
        j["eff_poly"] = TIME_F(eff_poly);
        j["Pout"] = TIME_F(Pout);
    }

#ifdef LIBTHERMO_USE_XTENSOR
    template <class G, class D = typename std::remove_reference_t<G>::value_t>
    void benchmark_vector(G&& gas, std::size_t size, unsigned long repeat, unsigned number)
    {
        typename xt::xtensor<D, 1>::shape_type shape = { size };
        xt::xtensor<D, 1> res(shape, 1.);

        xt::xtensor<D, 1> t1(shape, 273.15);
        xt::xtensor<D, 1> t2(shape, 350.);

        xt::xtensor<D, 1> p1(shape, 101325.);
        xt::xtensor<D, 1> p2(shape, 500000.);

        xt::xtensor<D, 1> eff(shape, 0.82);

        auto bench_cp = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.cp(t1);
            };
        };

        auto bench_gamma = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.gamma(t1);
            };
        };

        auto bench_h = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.h(t1);
            };
        };

        auto bench_phi = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.phi(t1);
            };
        };
        /*
            auto bench_dPhi = [&]() -> void {
                for (long i = 0; i < repeat; i++)
                {
                    xt::noalias(res) = gas.dphi(t1, t2);
                };
            };
        */
        /*
                auto bench_r = [&](const long &ntimes) -> void
                    { for(long i=0; i<ntimes; i++) { res = gas.template
           r<xt::xtensor<double, 1> >(); }; };
        */
        auto bench_pr = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.pr(t1, t2, eff);
            };
        };

        auto bench_eff_poly = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.eff_poly(p1, t1, p2, t2);
            };
        };

        auto bench_pout = [&]() -> void
        {
            for (long i = 0; i < repeat; i++)
            {
                xt::noalias(res) = gas.pr(t1, t2, gas.eff_poly(p1, t1, p2, t2)) * p2;
            };
        };

        nlohmann::json j;
        j["N"] = size;
        j["cp"] = timeit(bench_cp, size, repeat, number);
        j["gamma"] = timeit(bench_gamma, size, repeat, number);
        j["h"] = timeit(bench_h, size, repeat, number);
        j["phi"] = timeit(bench_phi, size, repeat, number);
        j["pr"] = timeit(bench_pr, size, repeat, number);
        j["eff_poly"] = timeit(bench_eff_poly, size, repeat, number);
        j["Pout"] = timeit(bench_pout, size, repeat, number);

        std::cout << j.dump() << std::endl;

        // fs::path results_f = "benchs/results/benchs_results.json";
        // nlohmann::json existing_j;
        // if (fs::exists(results_f))
        // {
        //     std::ifstream current_file(results_f);
        //     current_file >> existing_j;
        // }

        // nlohmann::json results;
        // results[gas.name()] = merge_results(existing_j[gas.name()], j);

        // std::ofstream out_file(results_f);
        // out_file << results.dump(2);
        // out_file.close();
    }
#endif
}  // namespace thermo

int
main()
{
    using namespace thermo;
    double Tref = 288.15;
    {
        IdealGas gas(287.05287, 1004.685045);

        std::cout << "\n"
                  << "Reference values" << std::endl;
        std::cout << "cp -> " << gas.cp(Tref) << std::endl;
        std::cout << "gamma -> " << gas.gamma(Tref) << std::endl;
        std::cout << "h -> " << gas.h(Tref) << std::endl;
        std::cout << "phi -> " << gas.phi(Tref) << std::endl;
        std::cout << "pr -> " << gas.pr(Tref, 450., 0.82) << std::endl;
        std::cout << "t_f_h -> " << gas.t_f_h(gas.h(Tref)) << std::endl;
        std::cout << "t_f_phi -> " << gas.t_f_phi(gas.phi(Tref)) << std::endl;
        std::cout << "t_f_pr -> " << gas.t_f_pr(gas.pr(Tref, 450., 0.82), Tref, 0.82) << std::endl;

        std::cout << "\nSingle value tests" << std::endl;
        // benchmark_single_value(gas, 1000000, 50);
        /*
                        std::cout << "\nLoop tests" << std::endl;
                        benchmark_loop(gas, 1000000, 500);

                #ifdef LIBTHERMO_USE_XTENSOR
                        std::cout << "\nVector tests" << std::endl;
                        benchmark_vector(gas, 1000000, 5000);
                #endif
                */
    }

    {
        PolyGas gas(PolyGasProps<double, 8>({ 0., 0., 0., 0., 0., 0., 0., 1004.4 }, 0., 287.05287));

        std::cout << "\n"
                  << "Reference values" << std::endl;
        std::cout << "cp -> " << gas.cp(Tref) << std::endl;
        std::cout << "gamma -> " << gas.gamma(Tref) << std::endl;
        std::cout << "h -> " << gas.h(Tref) << std::endl;
        std::cout << "phi -> " << gas.phi(Tref) << std::endl;
        std::cout << "pr -> " << gas.pr(Tref, 450., 0.82) << std::endl;
        std::cout << "t_f_h -> " << gas.t_f_h(gas.h(Tref), 1e-8) << std::endl;
        std::cout << "t_f_phi -> " << gas.t_f_phi(gas.phi(Tref), 1e-8) << std::endl;
        std::cout << "t_f_pr -> " << gas.t_f_pr(gas.pr(Tref, 450., 0.82), Tref, 0.82, 1e-8) << std::endl;

        std::cout << "\nSingle value tests" << std::endl;
        benchmark_single_value(gas, 100000, 50);

        // std::cout << "\nLoop tests" << std::endl;
        // benchmark_loop(gas, 10, 10000000);
        // benchmark_loop(gas, 100, 1000000);

#ifdef LIBTHERMO_USE_XTENSOR
        std::cout << "\nVector tests 1M" << std::endl;
        benchmark_vector(gas, 1000000, 10, 10);

        PolyGas gas_float(PolyGasProps<float, 8>({ 0., 0., 0., 0., 0., 0., 0., 1004.4 }, 0., 287.05287));
        std::cout << "\nVector tests 1M float" << std::endl;
        benchmark_vector(gas_float, 1000000, 10, 10);
#endif
    }
    // std::cout << XSIMD_X86_INSTR_SET << std::endl;

    return 0;
}
