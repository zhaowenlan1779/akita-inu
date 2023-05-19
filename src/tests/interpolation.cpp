// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include "core/interpolation.h"

TEST_CASE("univariate interpolation", "[DEPIR]") {
    UnivariateInterpolator interp{19};

    std::array<UnivariateInterpolator::dInt, 19> ys;
    for (std::size_t i = 0; i < ys.size(); ++i) {
        if (i % 3 == 0) {
            ys[i] = UnivariateInterpolator::dInt{static_cast<int64_t>(i + 3), 19};
        } else if (i % 3 == 1) {
            ys[i] = UnivariateInterpolator::dInt{static_cast<int64_t>(i * i + 2), 19};
        } else {
            ys[i] = UnivariateInterpolator::dInt{static_cast<int64_t>(2 * i + 1), 19};
        }
    }
    const auto& poly = interp.Interpolate(ys);

    // Evaluate poly
    const auto Evaluate = [&poly](UnivariateInterpolator::dInt x) {
        UnivariateInterpolator::dInt cur{1, 19};
        UnivariateInterpolator::dInt result{0, 19};
        for (std::size_t i = 0; i < poly.coeffs.size(); ++i) {
            result += poly.coeffs[i] * cur;
            cur *= x;
        }
        return result;
    };
    for (std::size_t i = 0; i < 19; ++i) {
        REQUIRE(Evaluate(UnivariateInterpolator::dInt{static_cast<int64_t>(i), 19}) == ys[i]);
    }
}
