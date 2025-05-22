#include "libthermo/poly_gas.hpp"

#include <iostream>
#include <chrono>
#include <array>
#include <stdexcept>


using namespace thermo;

static constexpr std::array<double, 8> air_cp_coeffs{
    2.35822e-20, -1.79613e-16, 4.70972e-13, -3.3094e-10, -6.27984e-07, 0.00123785, -0.521742, 1068.43,
};
static constexpr std::array<double, 8> fuel_cp_coeffs{ 0.0, 0.0, 0.0, 0.0, 0.0, -0.000527394, 2.3021158, 1285.0007 };
static constexpr std::array<double, 8> water_cp_coeffs{ 0.0,         0.0,        -3.86134e-14, 3.345e-10,
                                                        -1.14179e-6, 0.00173478, -0.484304,    1884.57 };

static constexpr std::array<double, 8> h_scalars{ 1. / 8., 1. / 7., 1. / 6., 1. / 5., 1. / 4., 1. / 3., 1. / 2., 1. };
static constexpr std::array<double, 8> dc_scalars{ 7., 6., 5., 4., 3., 2., 1., 1. };
static constexpr std::array<double, 8> fi_scalars{ 1. / 7., 1. / 6., 1. / 5., 1. / 4., 1. / 3., 1. / 2., 1., 1. };


template <class C, int S>
struct TP1906Props : PolyIdealGasMultiProps<C, S>
{
    using result_type = typename PolyIdealGasMultiProps<C, S>::result_type;

    double m_far = 0., m_war = 0.;

    TP1906Props() = default;

    TP1906Props(result_type FARB, result_type WAR)
    {
        update(FARB, WAR);
    }

    void update(result_type FARB, result_type WAR)
    {
        if (FARB.size() != WAR.size())
        {
            std::cout << "Size error" << std::endl;
            throw std::runtime_error("Invalid");
        }

        this->cp_coeffs.resize(FARB.size() * 8);
        this->dcp_coeffs.resize(FARB.size() * 8);  // Allocate 8 instead of 7 slots to allow SIMD
        this->h_coeffs.resize(FARB.size() * 9);
        this->phi_coeffs.resize(FARB.size() * 8);

        this->r.resize(FARB.size());
        this->phi_log.resize(FARB.size());

        using avx2_type = xsimd::batch<double, xsimd::avx2>;
        std::size_t inc = avx2_type::size;

        std::size_t size = FARB.size();

        for (std::size_t b = 0; b < size; b += inc)
        {
            avx2_type farb_vec = avx2_type::load_unaligned(&FARB[b]);
            avx2_type war_vec = avx2_type::load_unaligned(&WAR[b]);
            avx2_type total = 1. / (1. + farb_vec + war_vec);

            for (std::size_t i = 0; i < 8; ++i)
            {
                avx2_type air_cp_coeffs_vec = air_cp_coeffs[i];
                avx2_type fuel_cp_coeffs_vec = fuel_cp_coeffs[i];
                avx2_type water_cp_coeffs_vec = water_cp_coeffs[i];

                avx2_type h_scalars_vec = h_scalars[i];
                avx2_type dc_scalars_vec = dc_scalars[i];
                avx2_type fi_scalars_vec = fi_scalars[i];

                avx2_type cp_mix_vec = xsimd::fma(fuel_cp_coeffs_vec,
                                                  farb_vec,
                                                  xsimd::fma(water_cp_coeffs_vec, war_vec, air_cp_coeffs_vec))
                                       * total;

                cp_mix_vec.store_unaligned(&this->cp_coeffs[i * size + b]);

                avx2_type dc_mix_vec = cp_mix_vec * dc_scalars_vec;
                dc_mix_vec.store_unaligned(&this->dcp_coeffs[i * size + b]);

                avx2_type fi_mix_vec = cp_mix_vec * fi_scalars_vec;
                fi_mix_vec.store_unaligned(&this->phi_coeffs[i * size + b]);

                avx2_type h_mix_vec = cp_mix_vec * h_scalars_vec;
                h_mix_vec.store_unaligned(&this->h_coeffs[i * size + b]);
            }

            avx2_type r_vec = (287.05287 + farb_vec * 287.05287 + war_vec * 461.522) * total;
            r_vec.store_unaligned(&this->r[b]);

            avx2_type cp7_vec = avx2_type::load_unaligned(&this->cp_coeffs[7 * size + b]);
            cp7_vec.store_unaligned(&this->phi_log[b]);

            avx2_type h0 = -42785.8 * farb_vec * total;
            h0.store_unaligned(&this->h_coeffs[8 * size + b]);
        }
    }
};


int
main()
{
    constexpr int size = 1024;
    constexpr int repeat = 100000;
    auto t1 = std::chrono::high_resolution_clock::now();

    std::vector<double> res(size);
    std::vector<double> farb_ref(size);
    std::vector<double> war_ref(size);

    TP1906Props<std::vector<double>, size> foo(farb_ref, war_ref);

    MultiPolyIdealGas g1(TP1906Props<std::vector<double>, size>(farb_ref, war_ref));

    auto& props1 = g1.properties();

    for (auto i = 0; i < repeat; ++i)
    {
        props1.update(farb_ref, war_ref);
        g1.h(288.15, res);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto t = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
             / static_cast<double>(repeat);
    std::cout << t << std::endl;

    g1.properties().update(farb_ref, war_ref);
    // std::cout << "libthermo r: " << g1.r() << " J/kg\n";
    // std::cout << "libthermo cp: " << g1.cp(288.15) << std::endl;
    g1.h(288.15, res);
    std::cout << "libthermo h: " << res[0] << std::endl;
    std::cout << "libthermo h: " << res[1] << std::endl;
    std::cout << "libthermo h: " << res[2] << std::endl;
    std::cout << "libthermo h: " << res[3] << std::endl;
    std::cout << "libthermo h: " << res[4] << std::endl;
    std::cout << "libthermo h: " << res[10] << std::endl;
    std::cout << "libthermo h: " << res[11] << std::endl;

    return 0;
}