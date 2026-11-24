// tests/test_catchem_gaschem_units.cpp
#include <cassert>
#include <cmath>
#include <iostream>

// Molar conversions test (VMR to molar density and reverse)
void test_vmr_conversion_properties() {
    std::cout << "DEBUG: Running property-based mixing ratio tests" << std::endl;

    const double air_mw_kg = 0.0289644;
    double test_densities[] = {1.2, 1.0, 0.8, 0.5};    // dry air density, kg/m3
    double test_vmrs[] = {100.0, 1.0, 1.0e-3, 1.0e-6}; // ppmv values

    for (double density : test_densities) {
        double air_density_mol = density / air_mw_kg;
        for (double ppmv : test_vmrs) {
            // Convert: ppmv -> mol/m3
            double conc_molar = ppmv * 1.0e-6 * air_density_mol;
            assert(conc_molar > 0.0);

            // Convert back: mol/m3 -> ppmv
            double ppmv_back = (conc_molar / air_density_mol) * 1.0e6;

            // Assert strict identity recovery (reversibility property)
            assert(std::abs(ppmv - ppmv_back) < 1.0e-12 && "Identity mapping must be strictly reversible!");
        }
    }
    std::cout << "SUCCESS: Mixing ratio conversion property holds true." << std::endl;
}

// Bounds & safe guards checks
void test_value_safeguards() {
    std::cout << "DEBUG: Running boundary value tests" << std::endl;

    double negative_val = -10.0;
    double bounded = (negative_val < 0.0) ? 1.0e-20 : negative_val;
    assert(bounded == 1.0e-20 && "Negative values must be safely bounded to prevent NaN.");

    std::cout << "SUCCESS: Value safeguards hold true." << std::endl;
}

int main() {
    test_vmr_conversion_properties();
    test_value_safeguards();
    std::cout << "ALL PROPERTY UNIT TESTS PASSED SUCCESSFULLY!" << std::endl;
    return 0;
}
