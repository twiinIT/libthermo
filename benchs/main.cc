#include "libthermo/poly_gas.hpp"

#include <memory>
#include <iostream>
#include <chrono>
#include <array>


using namespace thermo;

static constexpr std::array<double, 8> air_cp_coeffs{
    2.35822e-20, -1.79613e-16, 4.70972e-13, -3.3094e-10, -6.27984e-07, 0.00123785, -0.521742, 1068.43,
};
static constexpr std::array<double, 8> fuel_cp_coeffs{ 0.0, 0.0, 0.0, 0.0, 0.0, -0.000527394, 2.3021158, 1285.0007 };
static constexpr std::array<double, 8> water_cp_coeffs{ 0.0,         0.0,        -3.86134e-14, 3.345e-10,
                                                        -1.14179e-6, 0.00173478, -0.484304,    1884.57 };

// static constexpr std::array<double, 3> h0_coeffs{ 0.0, -42785.8, 0.0 };
// static constexpr std::array<double, 3> r_coeffs{ 287.05287, 287.05287, 461.522 };

static constexpr std::array<double, 8> h_scalars{ 1. / 8., 1. / 7., 1. / 6., 1. / 5., 1. / 4., 1. / 3., 1. / 2., 1. };
static constexpr std::array<double, 8> dc_scalars{ 7., 6., 5., 4., 3., 2., 1., 1. };
static constexpr std::array<double, 8> fi_scalars{ 1. / 7., 1. / 6., 1. / 5., 1. / 4., 1. / 3., 1. / 2., 1., 1. };


template <class T>
struct PolyGasDynamicProps
{
    using arr_t = std::vector<T>;
    using value_t = T;
    static constexpr int size = -1;

    struct PolyGasDynamicProps()
    {
        cp_coeffs.resize(8);
        h_coeffs.resize(9);
    }

    arr_t cp_coeffs;
    arr_t dcp_coeffs;
    arr_t h_coeffs;
    arr_t phi_coeffs;
    T r, phi_log;

    void update(double FARB, double WAR)
    {
        double total = 1. / (1. + FARB + WAR);

        for (auto i = 0; i < cp_coeffs.size(); ++i)
            cp_coeffs[i] = (air_cp_coeffs[i] + FARB * fuel_cp_coeffs[i] + WAR * water_cp_coeffs[i]) * total;

        for (auto i = 0; i < cp_coeffs.size(); ++i)
            h_coeffs[i] = cp_coeffs[i] * h_scalars[i];
        h_coeffs[8] = FARB * -42785.8 * total;
        r = 287.05287;
    }
};

template <class T>
struct TP1906Props : PolyIdealGasProps<double, 8>
{
    double m_far = 0., m_war = 0.;
    TP1906Props()
    {
        update(0, 0);
    }

    TP1906Props(double FARB, double WAR)
    {
        update(FARB, WAR);
    }

    void update(double FARB, double WAR)
    {
        // if (FARB + WAR == m_far + m_war)
        //     return;

        m_far = FARB;
        m_war = WAR;
        double total = 1. / (1. + FARB + WAR);

        r = (287.05287 + FARB * 287.05287 + WAR * 461.522) * total;

        using avx2_type = xsimd::batch<double, xsimd::avx2>;
        std::size_t inc = avx2_type::size;

        avx2_type farb_vec = FARB;
        avx2_type war_vec = WAR;

        for (std::size_t i = 0; i < 8; i += inc)
        {
            avx2_type air_cp_coeffs_vec = avx2_type::load_unaligned(&air_cp_coeffs[i]);
            avx2_type fuel_cp_coeffs_vec = avx2_type::load_unaligned(&fuel_cp_coeffs[i]);
            avx2_type water_cp_coeffs_vec = avx2_type::load_unaligned(&water_cp_coeffs[i]);

            avx2_type cp_mix_vec
                = xsimd::fma(fuel_cp_coeffs_vec, farb_vec, xsimd::fma(water_cp_coeffs_vec, war_vec, air_cp_coeffs_vec))
                  * total;

            cp_mix_vec.store_unaligned(&cp_coeffs[i]);

            avx2_type dc_scalars_vec = avx2_type::load_unaligned(&dc_scalars[i]);
            avx2_type dc_mix_vec = cp_mix_vec * dc_scalars_vec;
            dc_mix_vec.store_unaligned(&dcp_coeffs[i]);

            avx2_type fi_scalars_vec = avx2_type::load_unaligned(&fi_scalars[i]);
            avx2_type fi_mix_vec = cp_mix_vec * fi_scalars_vec;
            fi_mix_vec.store_unaligned(&phi_coeffs[i]);

            avx2_type h_scalars_vec = avx2_type::load_unaligned(&h_scalars[i]);
            avx2_type h_mix_vec = cp_mix_vec * h_scalars_vec;
            h_mix_vec.store_unaligned(&h_coeffs[i]);
            h_coeffs[8] = -42785.8 * FARB * total;

            phi_log = cp_coeffs[7];
        }
    }
};


int
main()
{
    constexpr int repeat = 1000000;
    volatile double res;
    PolyIdealGas<TP1906Props<double>> g{};

    auto t1 = std::chrono::high_resolution_clock::now();

    for (auto i = 0; i < repeat; ++i)
    {
        g.properties().update(0., 0.);
        res = g.h(288.15);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto t = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
             / static_cast<double>(repeat);
    std::cout << t << std::endl;

    constexpr int repeat2 = 100000;
    std::vector<PolyIdealGas<TP1906Props<double>>> gases(1024);
    t1 = std::chrono::high_resolution_clock::now();
    for (auto i = 0; i < repeat2; ++i)
    {
        for (auto& g : gases)
        {
            g.properties().update(0., 0.);
            res = g.h(288.15);
        }
    }
    t2 = std::chrono::high_resolution_clock::now();
    t = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
        / static_cast<double>(repeat2 * gases.size());
    std::cout << t << std::endl;

    // g1.properties().update(0.012, 0.0);
    // std::cout << "libthermo r: " << g1.r() << " J/kg\n";
    // std::cout << "libthermo cp: " << g1.cp(288.15) << std::endl;
    // std::cout << "libthermo h: " << g1.h(288.15) << std::endl;
    // std::cout << g1.pr(288.15, 400., 0.8) << std::endl;

    // g2.properties().update(0.012, 0.005);
    // std::cout << g2.cp(288.15) << std::endl;
    // std::cout << g2.h(288.15) << std::endl;

    return 0;
}