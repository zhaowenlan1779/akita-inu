// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <array>
#include <cmath>
#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include "core/ASHE.h"
#include "core/kissfft.hpp"

TEST_CASE("FFT", "[ASHE]") {
    SECTION("Recovery of coefficients") {
        // Quick verification that my maths is not wrong
        static constexpr std::size_t N = 5;

        std::array<std::complex<double>, 2 * N> A{};
        for (std::size_t i = 0; i < N; ++i) {
            A[i] = {static_cast<double>(i), static_cast<double>(i * i + 1)};
        }

        std::array<std::complex<double>, 2 * N> B{};

        kissfft<double> fft(2 * N, false);
        fft.transform(A.data(), B.data());

        for (std::size_t i = 0; i < 2 * N; i += 2) {
            B[i] = 0;
        }

        kissfft<double> inv(2 * N, true);
        inv.transform(B.data(), A.data());

        static constexpr double tol = 1e-6;
        for (std::size_t i = 0; i < N; ++i) {
            REQUIRE(std::abs(A[i].real() / N - i) < tol);
            REQUIRE(std::abs(A[i].imag() / N - (i * i + 1)) < tol);
        }
    }

    SECTION("Coefficients real when generated") {
        static constexpr std::size_t N = 10;
        std::array<std::complex<double>, 2 * N> A{};
        for (int i = 0; i < static_cast<int>(N / 2); ++i) {
            A[2 * i + 1] = {static_cast<double>(i), static_cast<double>(i * i + 1)};
            A[2 * (N - 1 - i) + 1] = {static_cast<double>(i), static_cast<double>(-i * i - 1)};
        }

        kissfft<double> fft(2 * N, true);

        std::array<std::complex<double>, 2 * N> B{};
        fft.transform(A.data(), B.data());

        static constexpr double tol = 1e-6;
        for (std::size_t i = 0; i < N; ++i) {
            REQUIRE(std::abs(B[i].imag()) < tol);
        }
    }
}

TEST_CASE("ASHE Correctness", "[ASHE]") {
    std::cout << "Performing small ASHE test, this may take a bit..." << std::endl;
    static constexpr std::size_t N = 10;

    static constexpr int64_t d = 91;
    ASHE ashe{64, d, 6, N};

    const auto s = ashe.Gen();

    static constexpr std::size_t m = 30;
    std::array<dInt, m> u;
    for (std::size_t i = 0; i < m; ++i) {
        u[i] = dInt{static_cast<int64_t>(i * i - 3), d};
    }

    std::array<ASHE::RElement, m> ct;
    for (std::size_t i = 0; i < m; ++i) {
        ct[i] = ashe.Enc(s, u[i]);
    }

    // Define the polynomial with 10 terms:
    // i-th term = (i * i + 2)X_i^2 X_{2i} X_{i + 20} X_{i^2 mod m}
    ASHE::RElement result{ashe.n, ashe.D};
    for (std::size_t i = 0; i < N; ++i) {
        result += ashe.Lift(dInt{static_cast<int64_t>(i * i + 2), d}) * ct[i] * ct[i] * ct[2 * i] *
                  ct[i + 20] * ct[(i * i) % m];
    }
    const auto decrypted = ashe.Dec(s, result);
    dInt answer{0, d};
    for (std::size_t i = 0; i < N; ++i) {
        answer += dInt{static_cast<int64_t>(i * i + 2), d} * u[i] * u[i] * u[2 * i] * u[i + 20] *
                  u[(i * i) % m];
    }
    REQUIRE(decrypted == answer);
}
