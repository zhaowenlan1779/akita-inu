// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <array>
#include <cmath>
#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include <kiss_fft.h>
#include "core/ASHE.h"

TEST_CASE("FFT", "[ASHE]") {
    SECTION("Recovery of coefficients") {
        // Quick verification that my maths is not wrong
        static constexpr std::size_t N = 5;

        std::array<kiss_fft_cpx, 2 * N> A;
        for (std::size_t i = 0; i < N; ++i) {
            A[i].r = i;
            A[i].i = i * i + 1;
        }

        std::array<kiss_fft_cpx, 2 * N> B;

        auto* cfg = kiss_fft_alloc(2 * N, false, 0, 0);
        kiss_fft(cfg, A.data(), B.data());
        kiss_fft_free(cfg);

        for (std::size_t i = 0; i < 2 * N; i += 2) {
            B[i].r = B[i].i = 0;
        }

        cfg = kiss_fft_alloc(2 * N, true, 0, 0);
        kiss_fft(cfg, B.data(), A.data());
        kiss_fft_free(cfg);

        static constexpr double tol = 1e-6;
        for (std::size_t i = 0; i < N; ++i) {
            REQUIRE(std::abs(A[i].r / N - i) < tol);
            REQUIRE(std::abs(A[i].i / N - (i * i + 1)) < tol);
        }
    }

    SECTION("Coefficients real when generated") {
        static constexpr std::size_t N = 10;
        std::array<kiss_fft_cpx, 2 * N> A{};
        for (int i = 0; i < static_cast<int>(N / 2); ++i) {
            A[2 * i + 1].r = i;
            A[2 * i + 1].i = i * i + 1;
            A[2 * (N - 1 - i) + 1].r = i;
            A[2 * (N - 1 - i) + 1].i = -(i * i + 1);
        }

        auto* cfg = kiss_fft_alloc(2 * N, true, 0, 0);
        kiss_fft(cfg, A.data(), A.data());
        kiss_fft_free(cfg);

        static constexpr double tol = 1e-6;
        for (std::size_t i = 0; i < N; ++i) {
            REQUIRE(std::abs(A[i].i) < tol);
        }
    }
}

TEST_CASE("ASHE Correctness", "[ASHE]") {
    std::cout << "Performing small ASHE test, this may take a bit..." << std::endl;

    auto d = std::make_shared<BigInt>(91);

    static constexpr std::size_t N = 10;

    ASHE ashe{64, d, 6, N};

    const auto s = ashe.Gen();

    static constexpr std::size_t m = 30;
    std::array<ASHE::dInt, m> u;
    for (std::size_t i = 0; i < m; ++i) {
        u[i] = ASHE::dInt{i * i - 3, d};
    }

    std::array<ASHE::RElement, m> ct;
    for (std::size_t i = 0; i < m; ++i) {
        ct[i] = ashe.Enc(s, u[i]);
    }

    // Define the polynomial with 10 terms:
    // i-th term = (i * i + 2)X_i^2 X_{2i} X_{i + 20} X_{i^2 mod m}
    ASHE::RElement result{ashe.n, ashe.D};
    for (std::size_t i = 0; i < N; ++i) {
        result += ashe.Lift(ASHE::dInt{i * i + 2, d}) * ct[i] * ct[i] * ct[2 * i] * ct[i + 20] *
                  ct[(i * i) % m];
    }
    const auto decrypted = ashe.Dec(s, result);
    ASHE::dInt answer{0, d};
    for (std::size_t i = 0; i < N; ++i) {
        answer += ASHE::dInt{i * i + 2, d} * u[i] * u[i] * u[2 * i] * u[i + 20] * u[(i * i) % m];
    }
    REQUIRE(decrypted == answer);
}
