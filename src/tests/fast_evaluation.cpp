// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include "core/fast_evaluation.h"

TEST_CASE("multidimensional FFT", "[DEPIR]") {
    MultiPoly poly{23, 3};
    for (int64_t i = 0; i < 23; ++i) {
        for (int64_t j = 0; j < 23; ++j) {
            for (int64_t k = 0; k < 23; ++k) {
                const std::size_t idx = i * 23 * 23 + j * 23 + k;
                if (idx % 3 == 0) {
                    poly.coeffs[idx] = dInt{i + j + k + 3, 23};
                } else if (idx % 3 == 1) {
                    poly.coeffs[idx] = dInt{i * j + 2, 23};
                } else {
                    poly.coeffs[idx] = dInt{j * k + 1, 23};
                }
            }
        }
    }

    Evaluation::MultidimensionalFFT fft{23, 3};
    const auto& result = fft.EvaluateAll(poly);
    const auto Evaluate = [&poly](dInt x1, dInt x2, dInt x3) {
        dInt result{0, 23};
        dInt c1{1, 23};
        for (std::size_t i = 0; i < 23; ++i) { // deg x1
            dInt c2{1, 23};
            for (std::size_t j = 0; j < 23; ++j) { // deg x2
                dInt c3{1, 23};
                for (std::size_t k = 0; k < 23; ++k) { // deg x3
                    const auto idx = k * 23 * 23 + j * 23 + i;
                    result += poly.coeffs[idx] * c1 * c2 * c3;
                    c3 *= x3;
                }
                c2 *= x2;
            }
            c1 *= x1;
        }
        return result;
    };

    for (int64_t i = 0; i < 23; ++i) {
        for (int64_t j = 0; j < 23; ++j) {
            for (int64_t k = 0; k < 23; ++k) {
                const std::size_t idx = i * 23 * 23 + j * 23 + k;
                REQUIRE(Evaluate({i, 23}, {j, 23}, {k, 23}) == result[idx]);
            }
        }
    }
}
