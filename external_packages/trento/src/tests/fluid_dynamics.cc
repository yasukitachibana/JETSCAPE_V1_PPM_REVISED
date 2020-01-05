#include "../fluid_dynamics.h"
#include "gtest/gtest.h"

void test_not_in_range(EvolutionHistory hist, real tau, real x, real y, real eta) {
    try {
        hist.check_in_range(tau, x, y, eta);
    } catch (InvalidSpaceTimeRange & e) {
        auto estr = std::string(e.what());
        ASSERT_TRUE(estr.find("not in range") != estr.npos);
    }
}

// test EvolutionHistory Class
TEST(EvolutionHistoryTest, TEST_WRITE){
    auto hist = EvolutionHistory();
    hist.tau_min = 0.6;
    hist.dtau = 0.1;
    hist.x_min = -10;
    hist.y_min = -10;
    hist.eta_min = -10;
    hist.dx = 0.1;
    hist.dy = 0.1;
    hist.deta = 0.1;
    hist.ntau = 10;
    hist.nx = 200;
    hist.ny = 200;
    hist.neta = 200;
    hist.tau_eta_is_tz = false;

    EXPECT_EQ(hist.tau_max(), static_cast<real>(1.6));
    EXPECT_EQ(hist.x_max(), static_cast<real>(10.0));
    EXPECT_EQ(hist.y_max(), static_cast<real>(10.0));
    EXPECT_EQ(hist.eta_max(), static_cast<real>(10.0));

    // check range test
    test_not_in_range(hist, 0.5, 0.3, 0.3, 0.3);
    test_not_in_range(hist, 20.5, 0.3, 0.3, 0.3);
    test_not_in_range(hist, 10.5, 10.3, 0.3, 0.3);
    test_not_in_range(hist, 10.5, 10., 10.3, 0.3);
    test_not_in_range(hist, 10.5, 10., 10., 10.3);

    real const_ed = 0.8;
    for (int n=0; n != hist.ntau; n++)
        for (int i=0; i != hist.nx; i++)
            for (int j=0; j != hist.ny; j++)
                for (int k=0; k != hist.neta; k++) {
                    auto cell = FluidCellInfo();
                    cell.energy_density = const_ed;
                    hist.data.emplace_back(std::move(cell));
                }
    EXPECT_EQ(hist.get(0.8, 0.0, 0.0, 0.0).energy_density, const_ed);
}
