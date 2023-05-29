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

TEST_CASE("Simple evaluator", "[DEPIR]") {
    using Evaluation::pInt;

    MultivariatePolynomial<pInt> poly{2, 2};
    for (auto& it : poly.coeffs) {
        it = pInt{static_cast<int16_t>(std::rand() % 20), 20};
    }

    const auto Evaluate = [&poly](const pInt& x1, const pInt& x2) {
        pInt result{0, 20};
        pInt c1{1, 20};
        for (std::size_t i = 0; i < 2; ++i) { // deg x1
            pInt c2{1, 20};
            for (std::size_t j = 0; j < 2; ++j) { // deg x2
                const auto idx = j * 2 + i;
                result += poly.coeffs[idx] * c1 * c2;
                c2 *= x2;
            }
            c1 *= x1;
        }
        return result;
    };

    static const auto TempPath = std::filesystem::u8path("tmp-simple-evaluator/");
    {
        // Save
        Evaluation::SimpleFastEvaluator eval{TempPath, poly, 20};
        for (std::size_t i = 0; i < 20; ++i) {
            for (std::size_t j = 0; j < 20; ++j) {
                const pInt x{static_cast<int16_t>(i), 20};
                const pInt y{static_cast<int16_t>(j), 20};
                const std::array<pInt, 2> xs{{x, y}};
                REQUIRE(eval.Evaluate(xs) == Evaluate(x, y));
            }
        }
    }
    {
        // Load
        Evaluation::SimpleFastEvaluator eval{TempPath};
        for (std::size_t i = 0; i < 20; ++i) {
            for (std::size_t j = 0; j < 20; ++j) {
                const pInt x{static_cast<int16_t>(i), 20};
                const pInt y{static_cast<int16_t>(j), 20};
                const std::array<pInt, 2> xs{{x, y}};
                REQUIRE(eval.Evaluate(xs) == Evaluate(x, y));
            }
        }
    }

    std::filesystem::remove_all(TempPath);
}

TEST_CASE("Scalar evaluator", "[DEPIR]") {
    const auto r = std::make_shared<BigInt>(100);

    using Evaluation::rInt;

    MultivariatePolynomial<rInt> poly{2, 2};
    for (auto& it : poly.coeffs) {
        it = rInt{std::rand() % 100, r};
    }

    const auto Evaluate = [&r, &poly](const rInt& x1, const rInt& x2) {
        rInt result{0, r};
        rInt c1{1, r};
        for (std::size_t i = 0; i < 2; ++i) { // deg x1
            rInt c2{1, r};
            for (std::size_t j = 0; j < 2; ++j) { // deg x2
                const auto idx = j * 2 + i;
                result += poly.coeffs[idx] * c1 * c2;
                c2 *= x2;
            }
            c1 *= x1;
        }
        return result;
    };

    static const auto TempPath = std::filesystem::u8path("tmp-scalar-evaluator/");
    {
        // Save
        Evaluation::ScalarFastEvaluator eval{TempPath, poly, r};
        for (std::size_t i = 0; i < 100; ++i) {
            const rInt x{std::rand() % 100, r};
            const rInt y{std::rand() % 100, r};
            const std::array<rInt, 2> xs{{x, y}};
            REQUIRE(eval.Evaluate(xs) == Evaluate(x, y));
        }
    }
    {
        // Load
        Evaluation::ScalarFastEvaluator eval{TempPath};
        for (std::size_t i = 0; i < 100; ++i) {
            const rInt x{std::rand() % 100, r};
            const rInt y{std::rand() % 100, r};
            const std::array<rInt, 2> xs{{x, y}};
            REQUIRE(eval.Evaluate(xs) == Evaluate(x, y));
        }
    }

    std::filesystem::remove_all(TempPath);
}

TEST_CASE("Univariate evaluator", "[DEPIR]") {
    const auto q = std::make_shared<BigInt>(3);

    using Element = Evaluation::UnivariateFastEvaluator::Element;
    const auto GenerateElement = [&q] {
        Element val{2};
        val.coeffs[0] = qInt{std::rand() % 3, q};
        val.coeffs[1] = qInt{std::rand() % 3, q};
        return val;
    };

    MultivariatePolynomial<Element> poly{2, 2};
    for (auto& it : poly.coeffs) {
        it = GenerateElement();
    }

    const auto Evaluate = [&poly](const Element& x1, const Element& x2) {
        Element result{2};

        Element c1{2};
        c1.coeffs[0] = BigInt{1};
        for (std::size_t i = 0; i < 2; ++i) { // deg x1
            Element c2{2};
            c2.coeffs[0] = BigInt{1};
            for (std::size_t j = 0; j < 2; ++j) { // deg x2
                const auto idx = j * 2 + i;
                result += poly.coeffs[idx] * c1 * c2;
                c2 *= x2;
            }
            c1 *= x1;
        }
        return result;
    };

    static const auto TempPath = std::filesystem::u8path("tmp-univariate-evaluator/");
    {
        // Save
        Evaluation::UnivariateFastEvaluator eval{TempPath, poly, q, 2};
        for (std::size_t i = 0; i < 100; ++i) {
            const auto x1 = GenerateElement();
            const auto x2 = GenerateElement();
            const std::array<Element, 2> xs{{x1, x2}};
            REQUIRE(eval.Evaluate(xs) == Evaluate(x1, x2));
        }
    }
    {
        // Load
        Evaluation::UnivariateFastEvaluator eval{TempPath};
        for (std::size_t i = 0; i < 100; ++i) {
            const auto x1 = GenerateElement();
            const auto x2 = GenerateElement();
            const std::array<Element, 2> xs{{x1, x2}};
            REQUIRE(eval.Evaluate(xs) == Evaluate(x1, x2));
        }
    }

    std::filesystem::remove_all(TempPath);
}

TEST_CASE("Bivariate evaluator", "[DEPIR]") {
    const auto q = std::make_shared<BigInt>(2);

    using Element = Evaluation::BivariateFastEvaluator::Element;
    const auto GenerateElement = [&q] {
        Element val{2, 2};
        for (auto& it : val.coeffs) {
            it = qInt{std::rand() % 2, q};
        }
        return val;
    };

    MultivariatePolynomial<Element> poly{2, 1};
    for (auto& it : poly.coeffs) {
        it = GenerateElement();
    }

    const auto Evaluate = [&poly](const Element& x1) {
        Element result{2, 2};

        Element c1{2, 2};
        c1.coeffs[0] = BigInt{1};
        for (std::size_t i = 0; i < 2; ++i) { // deg x1
            result += poly.coeffs[i] * c1;
            c1 *= x1;
        }
        return result;
    };

    static const auto TempPath = std::filesystem::u8path("tmp-bivariate-evaluator/");
    {
        // Save
        Evaluation::BivariateFastEvaluator eval{TempPath, poly, q, 2, 2};
        for (std::size_t i = 0; i < 10; ++i) {
            const auto x1 = GenerateElement();
            const std::array<Element, 1> xs{{x1}};
            REQUIRE(eval.Evaluate(xs) == Evaluate(x1));
        }
    }
    {
        // Load
        Evaluation::BivariateFastEvaluator eval{TempPath};
        for (std::size_t i = 0; i < 10; ++i) {
            const auto x1 = GenerateElement();
            const std::array<Element, 1> xs{{x1}};
            REQUIRE(eval.Evaluate(xs) == Evaluate(x1));
        }
    }

    std::filesystem::remove_all(TempPath);
}
