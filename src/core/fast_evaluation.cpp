// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <algorithm>
#include "core/fast_evaluation.h"
#include "core/kissfft.hpp"

namespace Evaluation {

static std::vector<int64_t> GetDistinctPrimeFactors(int64_t s) {
    int64_t cur = s;
    int64_t p = 2;
    std::vector<int64_t> factors;
    while (cur > 1) {
        while (cur % p) {
            if (p == 2) {
                p = 3;
            } else {
                p += 2;
            }
            if (p * p > s) {
                p = cur;
            }
        }
        factors.emplace_back(p);
        while (cur % p == 0) {
            cur /= p;
        }
    }
    return factors;
}

static int64_t QuickPow(int64_t base, int64_t power, int64_t modulus) {
    int64_t result = 1;
    while (power > 0) {
        if (power & 1) {
            result *= base;
            result %= modulus;
        }
        power >>= 1;
        base *= base;
        base %= modulus;
    }
    return result;
}

MultidimensionalFFT::MultidimensionalFFT(int64_t d_, std::size_t m_)
    : d(d_), m(m_), fft{static_cast<std::size_t>(d - 1), false} {

    const auto s = d - 1;
    const auto& factors = GetDistinctPrimeFactors(s);

    // Find the smallest primitive root
    for (int64_t t = 2; t < d; ++t) {
        if (std::all_of(factors.begin(), factors.end(),
                        [this, t, s](int64_t p) { return QuickPow(t, s / p, d) != 1; })) {
            alpha = t;
            break;
        }
    }
    if (!alpha) {
        throw std::runtime_error("Failed to find primitive root");
    }

    // Initialize twiddles
    fft._twiddles.resize(d - 1);
    dInt cur{1, d};
    for (int64_t i = 0; i < d - 1; ++i) {
        fft._twiddles[i] = cur;
        cur *= alpha;
    }
}

MultidimensionalFFT::~MultidimensionalFFT() = default;

std::vector<dInt> MultidimensionalFFT::EvaluateAll(const MultiPoly& poly) const {
    return EvaluateAllImpl(std::span{poly.coeffs}, m);
}

void MultidimensionalFFT::EvaluateUnivariate(std::span<const dInt> coeffs,
                                             std::span<dInt> out) const {
    std::vector<dInt> results(d - 1);
    fft.transform(coeffs.data(), results.data());

    int64_t cur = 1;
    for (std::size_t i = 0; i < results.size(); ++i) {
        out[cur] = results[i] + coeffs.back();
        cur *= alpha;
        cur %= d;
    }
    out[0] = coeffs[0];
}

std::vector<dInt> MultidimensionalFFT::EvaluateAllImpl(std::span<const dInt> coeffs,
                                                       std::size_t stage) const {
    std::vector<dInt> out(coeffs.size());
    if (stage == 1) {
        EvaluateUnivariate(coeffs, out);
        return out;
    }
    std::vector<dInt> final_coeffs(coeffs.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(d); ++i) {
        auto evals =
            EvaluateAllImpl({coeffs.begin() + coeffs.size() / d * i, coeffs.size() / d}, stage - 1);
        for (std::size_t j = 0; j < evals.size(); ++j) {
            final_coeffs[j * d + i] = evals[j];
        }
    }

    for (std::size_t i = 0; i < coeffs.size() / d; ++i) {
        EvaluateUnivariate({final_coeffs.begin() + d * i, static_cast<std::size_t>(d)},
                           {out.begin() + d * i, static_cast<std::size_t>(d)});
    }
    return out;
}

} // namespace Evaluation
