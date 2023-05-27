// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include "core/fast_evaluation.h"

TEST_CASE("multidimensional FFT", "[DEPIR]") {
    MultivariatePolynomial<Evaluation::piInt> poly{23, 3};
    for (int16_t i = 0; i < 23; ++i) {
        for (int16_t j = 0; j < 23; ++j) {
            for (int16_t k = 0; k < 23; ++k) {
                const std::size_t idx = i * 23 * 23 + j * 23 + k;
                if (idx % 3 == 0) {
                    poly.coeffs[idx] =
                        Evaluation::piInt{static_cast<int16_t>(i + j + k + 3), int16_t{17}};
                } else if (idx % 3 == 1) {
                    poly.coeffs[idx] =
                        Evaluation::piInt{static_cast<int16_t>(i * j + 2), int16_t{17}};
                } else {
                    poly.coeffs[idx] =
                        Evaluation::piInt{static_cast<int16_t>(j * k + 1), int16_t{17}};
                }
            }
        }
    }

    Evaluation::MultidimensionalFFT fft{17, 3};
    const auto& result = fft.EvaluateAll(poly);
    const auto Evaluate = [&poly](Evaluation::piInt x1, Evaluation::piInt x2,
                                  Evaluation::piInt x3) {
        Evaluation::piInt result{0, 17};
        Evaluation::piInt c1{1, 17};
        for (std::size_t i = 0; i < 23; ++i) { // deg x1
            Evaluation::piInt c2{1, 17};
            for (std::size_t j = 0; j < 23; ++j) { // deg x2
                Evaluation::piInt c3{1, 17};
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

    for (int16_t i = 0; i < 17; ++i) {
        for (int16_t j = 0; j < 17; ++j) {
            for (int16_t k = 0; k < 17; ++k) {
                const std::size_t idx = i * 17 * 17 + j * 17 + k;
                REQUIRE(Evaluate({i, 17}, {j, 17}, {k, 17}) == result[idx]);
            }
        }
    }
}

TEST_CASE("Simple evaluator", "[.DEPIR]") {
    MultivariatePolynomial<Evaluation::pInt> poly{3, 3};
    for (int16_t i = 0; i < 3; ++i) {
        for (int16_t j = 0; j < 3; ++j) {
            for (int16_t k = 0; k < 3; ++k) {
                const std::size_t idx = i * 3 * 3 + j * 3 + k;
                if (idx % 4 == 0) {
                    poly.coeffs[idx] = Evaluation::pInt{
                        static_cast<int16_t>((i + 2) * (j + 1) * (k + 4) + 3), int16_t{137}};
                } else if (idx % 4 == 1) {
                    poly.coeffs[idx] = Evaluation::pInt{
                        static_cast<int16_t>((i + 1) * (j + 2) + 2 * k + 2), int16_t{137}};
                } else {
                    poly.coeffs[idx] = Evaluation::pInt{
                        static_cast<int16_t>((j + 4) * (k + 6) + i + 1), int16_t{137}};
                }
            }
        }
    }

    Evaluation::SimpleFastEvaluator eval{std::filesystem::u8path("out/"), poly, 137};
    const auto Evaluate = [&poly](Evaluation::pInt x1, Evaluation::pInt x2, Evaluation::pInt x3) {
        Evaluation::pInt result{0, 137};
        Evaluation::pInt c1{1, 137};
        for (std::size_t i = 0; i < 3; ++i) { // deg x1
            Evaluation::pInt c2{1, 137};
            for (std::size_t j = 0; j < 3; ++j) { // deg x2
                Evaluation::pInt c3{1, 137};
                for (std::size_t k = 0; k < 3; ++k) { // deg x3
                    const auto idx = k * 3 * 3 + j * 3 + i;
                    result += poly.coeffs[idx] * c1 * c2 * c3;
                    c3 *= x3;
                }
                c2 *= x2;
            }
            c1 *= x1;
        }
        return result;
    };

    for (std::size_t i = 0; i < 10000; ++i) {
        const Evaluation::pInt x{static_cast<int16_t>(std::rand() % 137), 137};
        const Evaluation::pInt y{static_cast<int16_t>(std::rand() % 137), 137};
        const Evaluation::pInt z{static_cast<int16_t>(std::rand() % 137), 137};
        const std::array<Evaluation::pInt, 3> xs{{x, y, z}};
        REQUIRE(eval.Evaluate(xs) == Evaluate(x, y, z));
    }
}
