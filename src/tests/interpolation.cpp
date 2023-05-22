// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <array>
#include <catch2/catch_test_macros.hpp>
#include "core/interpolation.h"

TEST_CASE("univariate interpolation", "[DEPIR]") {
    Interpolation::UnivariateInterpolator interp{19};

    std::array<dInt, 19> ys;
    for (std::size_t i = 0; i < ys.size(); ++i) {
        if (i % 3 == 0) {
            ys[i] = dInt{static_cast<int64_t>(i + 3), 19};
        } else if (i % 3 == 1) {
            ys[i] = dInt{static_cast<int64_t>(i * i + 2), 19};
        } else {
            ys[i] = dInt{static_cast<int64_t>(2 * i + 1), 19};
        }
    }
    const auto& poly = interp.Interpolate(ys);

    // Evaluate poly
    const auto Evaluate = [&poly](dInt x) {
        dInt cur{1, 19};
        dInt result{0, 19};
        for (std::size_t i = 0; i < poly.coeffs.size(); ++i) {
            result += poly.coeffs[i] * cur;
            cur *= x;
        }
        return result;
    };
    for (std::size_t i = 0; i < 19; ++i) {
        REQUIRE(Evaluate(dInt{static_cast<int64_t>(i), 19}) == ys[i]);
    }
    REQUIRE(poly.coeffs.size() <= 19);
}

TEST_CASE("multivariate interpolation", "[DEPIR]") {
    Interpolation::MultivariateInterpolator interp{19, 3};

    // Generate the Ys
    std::array<dInt, 19 * 19 * 19> ys;
    for (int64_t i = 0; i < 19; ++i) {
        for (int64_t j = 0; j < 19; ++j) {
            for (int64_t k = 0; k < 19; ++k) {
                const std::size_t idx = i * 19 * 19 + j * 19 + k;
                if (idx % 3 == 0) {
                    ys[idx] = dInt{i + j + k + 3, 19};
                } else if (idx % 3 == 1) {
                    ys[idx] = dInt{i * j + 2, 19};
                } else {
                    ys[idx] = dInt{j * k + 1, 19};
                }
            }
        }
    }
    const auto& poly = interp.Interpolate(ys);
    const auto Evaluate = [&poly](dInt x1, dInt x2, dInt x3) {
        dInt result{0, 19};
        dInt c1{1, 19};
        for (std::size_t i = 0; i < 19; ++i) { // deg x1
            dInt c2{1, 19};
            for (std::size_t j = 0; j < 19; ++j) { // deg x2
                dInt c3{1, 19};
                for (std::size_t k = 0; k < 19; ++k) { // deg x3
                    const auto idx = k * 19 * 19 + j * 19 + i;
                    result += poly.coeffs[idx] * c1 * c2 * c3;
                    c3 *= x3;
                }
                c2 *= x2;
            }
            c1 *= x1;
        }
        return result;
    };

    for (int64_t i = 0; i < 19; ++i) {
        for (int64_t j = 0; j < 19; ++j) {
            for (int64_t k = 0; k < 19; ++k) {
                const std::size_t idx = i * 19 * 19 + j * 19 + k;
                REQUIRE(Evaluate({i, 19}, {j, 19}, {k, 19}) == ys[idx]);
            }
        }
    }
}
