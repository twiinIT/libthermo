#include "libthermo/ideal_gas.hpp"
#include "libthermo/poly_gas.hpp"

#include "cantera/thermo/NasaPoly1.h"
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/thermo/Species.h>


#include <memory>
#include <iostream>
#include <chrono>
#include <array>


using namespace Cantera;
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
    TP1906Props(double FARB, double WAR)
    {
        update(FARB, WAR);
    }

    void update(double FARB, double WAR)
    {
        if (FARB + WAR == m_far + m_war)
            return;

        m_far = FARB;
        m_war = WAR;
        double total = 1. / (1. + FARB + WAR);

        r = (287.05287 + FARB * 287.05287 + WAR * 461.522) * total;

        using avx2_type = xsimd::batch<double, xsimd::avx2>;
        // using sse2_type = xsimd::batch<double, xsimd::sse2>;
        std::size_t inc = avx2_type::size;

        avx2_type farb_vec = FARB;
        avx2_type war_vec = WAR;

        for (std::size_t i = 4; i < 8; i -= inc)
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


class CustomPoly final
    : public SpeciesThermoInterpType
    , public PolyIdealGas<PolyIdealGasProps<double, 8>>
{
public:
    CustomPoly()
    // : m_coeff(7)
    // , m_coeff5_orig(0.0)
    {
    }

    //! Constructor with all input data
    /*!
     * @param tlow    Minimum temperature
     * @param thigh   Maximum temperature
     * @param pref    reference pressure (Pa).
     * @param coeffs  Vector of coefficients used to set the parameters for the
     *                standard state, in the order [a0,a1,a2,a3,a4,a5,a6]
     */
    CustomPoly(double tlow, double thigh, double pref, std::array<double, 8> coeffs, double h0)
        : SpeciesThermoInterpType(tlow, thigh, pref)
    // , m_coeff(coeffs, coeffs + 7)
    {
        m_props = PolyIdealGasProps<double, 8>(coeffs, h0, 287.058);
        // m_coeff5_orig = m_coeff[5];
    }


    //! Set array of 7 polynomial coefficients
    void setParameters(double farb, double war)
    {
    }

    // int reportType() const override
    // {
    //     return NASA1;
    // }

    // size_t temperaturePolySize() const override
    // {
    //     return 8;
    // }

    void updateTemperaturePoly(double T, double* T_poly) const override
    {
        T_poly[0] = T;
        // T_poly[1] = T * T;
        // T_poly[2] = T_poly[1] * T;
        // T_poly[3] = T_poly[2] * T;
        // T_poly[4] = 1.0 / T;
        // T_poly[5] = std::log(T);
    }

    /**
     * @copydoc SpeciesThermoInterpType::updateProperties
     *
     * Temperature Polynomial:
     *  tt[0] = t;
     *  tt[1] = t*t;
     *  tt[2] = m_t[1]*t;
     *  tt[3] = m_t[2]*t;
     *  tt[4] = 1.0/t;
     *  tt[5] = std::log(t);
     */
    void updateProperties(const double* tt, double* cp_R, double* h_RT, double* s_R) const override
    {
        // double ct0 = m_coeff[0];          // a0
        // double ct1 = m_coeff[1] * tt[0];  // a1 * T
        // double ct2 = m_coeff[2] * tt[1];  // a2 * T^2
        // double ct3 = m_coeff[3] * tt[2];  // a3 * T^3
        // double ct4 = m_coeff[4] * tt[3];  // a4 * T^4

        // double cp, h, s;
        // cp = ct0 + ct1 + ct2 + ct3 + ct4;
        // h = ct0 + 0.5 * ct1 + 1.0 / 3.0 * ct2 + 0.25 * ct3 + 0.2 * ct4 + m_coeff[5] * tt[4];  // last term is a5/T
        // s = ct0 * tt[5] + ct1 + 0.5 * ct2 + 1.0 / 3.0 * ct3 + 0.25 * ct4 + m_coeff[6];        // last term is a6

        // return the computed properties for this species
        *cp_R = cp(*tt) / 290.901;
        *h_RT = h(*tt) / 290.901 / *tt;
        // *s_R = s;
    }

    // void updatePropertiesTemp(const double temp, double* cp_R, double* h_RT, double* s_R) const override
    // {
    //     double tPoly[6];
    //     updateTemperaturePoly(temp, tPoly);
    //     updateProperties(tPoly, cp_R, h_RT, s_R);
    // }

    // void reportParameters(
    //     size_t& n, int& type, double& tlow, double& thigh, double& pref, double* const coeffs) const override
    // {
    //     n = 0;
    //     type = NASA1;
    //     tlow = m_lowT;
    //     thigh = m_highT;
    //     pref = m_Pref;
    //     std::copy(m_coeff.begin(), m_coeff.end(), coeffs);
    // }

    // void getParameters(AnyMap& thermo) const override
    // {
    //     // NasaPoly1 is only used as an embedded model within NasaPoly2, so all
    //     // that needs to be added here are the polynomial coefficients
    //     thermo["data"].asVector<vector<double>>().push_back(m_coeff);
    // }

    // double reportHf298(double* const h298 = nullptr) const override
    // {
    //     double tt[6];
    //     double temp = 298.15;
    //     updateTemperaturePoly(temp, tt);
    //     double ct0 = m_coeff[0];          // a0
    //     double ct1 = m_coeff[1] * tt[0];  // a1 * T
    //     double ct2 = m_coeff[2] * tt[1];  // a2 * T^2
    //     double ct3 = m_coeff[3] * tt[2];  // a3 * T^3
    //     double ct4 = m_coeff[4] * tt[3];  // a4 * T^4

    //     double h_RT = ct0 + 0.5 * ct1 + 1.0 / 3.0 * ct2 + 0.25 * ct3 + 0.2 * ct4 + m_coeff[5] * tt[4];  // last t

    //     double h = h_RT * GasConstant * temp;
    //     if (h298)
    //     {
    //         *h298 = h;
    //     }
    //     return h;
    // }

    // void modifyOneHf298(const size_t k, const double Hf298New) override
    // {
    //     double hcurr = reportHf298(0);
    //     double delH = Hf298New - hcurr;
    //     m_coeff[5] += (delH) / GasConstant;
    // }

    void resetHf298() override
    {
        // m_coeff[5] = m_coeff5_orig;
    }

    // protected:
    //     //! array of polynomial coefficients, stored in the order [a0, ..., a6]
    //     vector<double> m_coeff;

    //     double m_coeff5_orig;
};

int
main()
{
    double tmin = 200.0;
    double tmax = 3500.0;
    double pref = 101325.;

    // NASA polynomial coefficients from GRI-Mech or Cantera's data
    // For simplicity, single range NasaPoly2 used for each species
    // Format: a[7] = {a1, a2, a3, a4, a5, a6, a7}

    auto dry_air = std::make_shared<Species>(
        "DryAir", Composition({ { "N", 0.7808 * 2 }, { "O", 0.2095 * 2 }, { "C", 0.0004 } }));
    // std::vector<double> dry_air_coeffs{
    //     1068.43 / 290.901, -0.521742 / 290.901, 0.00123785 / 290.901, -6.27984e-07 / 290.901, 0., 0., 0.
    // };
    // dry_air->thermo = std::make_shared<NasPoly1>(tmin, tmax, pref, dry_air_coeffs.data());
    dry_air->thermo = std::make_shared<CustomPoly>(tmin, tmax, pref, air_cp_coeffs, 0.);

    auto fuel = std::make_shared<Species>("FuelB",
                                          Composition({ { "N", 0.7808 * 2 }, { "O", 0.2095 * 2 }, { "C", 0.0004 } }));
    // std::vector<double> fuel_coeffs{
    //     1285.0007 / 290.901, 2.3021158 / 290.901, -0.000527394 / 290.901, 0., 0., -42785.8 / 290.901, 0.,
    // };
    // fuel->thermo = std::make_shared<NasPoly1>(tmin, tmax, pref, fuel_coeffs.data());
    fuel->thermo = std::make_shared<CustomPoly>(tmin, tmax, pref, fuel_cp_coeffs, -42785.8);

    // Create IdealGasPhase and initialize
    IdealGasPhase cantera_gas;
    cantera_gas.addSpecies(dry_air);
    cantera_gas.addSpecies(fuel);
    cantera_gas.initThermo();

    // Set temperature, pressure, and composition of dry air
    cantera_gas.setState_TP(288.15, OneAtm);

    constexpr int repeat = 10000000;
    auto t1 = std::chrono::high_resolution_clock::now();

    volatile double res;

    PolyIdealGas g1(TP1906Props<double>(0., 0.));
    PolyIdealGas<PolyGasDynamicProps<double>> g2;

    auto& props1 = g1.properties();
    auto& props2 = g2.properties();
    std::array<double, 2> composition{ 1, 0 };

    for (auto i = 0; i < repeat; ++i)
    {
        // composition[1] = i * 1e-10;
        // cantera_gas.setTemperature(288.15 + i * 1e-10);
        // cantera_gas.setMassFractions(composition.data());
        // // cantera_gas.setState_TPX(288.15, 101325., composition.data());
        // res = cantera_gas.enthalpy_mass();

        props1.update(i * 1e-10, 0.);
        res = g1.h(288.15);
        // res = g1.pr(288.15, 400., 0.8);

        // props2.update(i * 1e-12, 0.);
        // res = g2.h(288.15);
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    auto t = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count())
             / static_cast<double>(repeat);
    std::cout << t << std::endl;


    // Evaluate and print enthalpy
    composition = { 1, 0.012 };
    cantera_gas.setMassFractions(composition.data());
    cantera_gas.setTemperature(288.15);
    std::cout << "cantera r: " << (GasConstant / cantera_gas.meanMolecularWeight()) << " J/kg\n";
    std::cout << "cantera cp: " << cantera_gas.cp_mass() << " J/kg\n";
    std::cout << "cantera h: " << cantera_gas.enthalpy_mass() << " J/kg\n";

    g1.properties().update(0.012, 0.0);
    std::cout << "libthermo r: " << g1.r() << " J/kg\n";
    std::cout << "libthermo cp: " << g1.cp(288.15) << std::endl;
    std::cout << "libthermo h: " << g1.h(288.15) << std::endl;
    // std::cout << g1.pr(288.15, 400., 0.8) << std::endl;

    // g2.properties().update(0.012, 0.005);
    // std::cout << g2.cp(288.15) << std::endl;
    // std::cout << g2.h(288.15) << std::endl;

    return 0;
}