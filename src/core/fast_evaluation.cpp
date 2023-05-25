// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <algorithm>
#include <fstream>
#include <type_traits>
#include <boost/integer/extended_euclidean.hpp>
#include "core/fast_evaluation.h"
#include "core/kissfft.hpp"

namespace Evaluation {

static std::vector<int16_t> GetDistinctPrimeFactors(int16_t s) {
    int16_t cur = s;
    int16_t p = 2;
    std::vector<int16_t> factors;
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

static int16_t QuickPow(int16_t base, int16_t power, int16_t modulus) {
    int16_t result = 1;
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

MultidimensionalFFT::MultidimensionalFFT(int16_t p_, std::size_t m_)
    : p(p_), m(m_), fft{static_cast<std::size_t>(p - 1), false} {

    if (p == 2) {
        return;
    }

    const auto s = p - 1;
    const auto& factors = GetDistinctPrimeFactors(s);

    // Find the smallest primitive root
    for (int16_t t = 2; t < p; ++t) {
        if (std::all_of(factors.begin(), factors.end(),
                        [this, t, s](int16_t q) { return QuickPow(t, s / q, p) != 1; })) {
            alpha = t;
            break;
        }
    }
    if (!alpha) {
        throw std::runtime_error("Failed to find primitive root");
    }

    // Initialize twiddles
    fft._twiddles.resize(p - 1);
    piInt cur{1, p};
    for (int16_t i = 0; i < p - 1; ++i) {
        fft._twiddles[i] = cur;
        cur *= alpha;
    }
}

MultidimensionalFFT::~MultidimensionalFFT() = default;

std::vector<piInt> MultidimensionalFFT::EvaluateAll(
    const MultivariatePolynomial<piInt>& poly) const {
    if (poly.d == static_cast<std::size_t>(p)) {
        return EvaluateAllImpl(std::span{poly.coeffs}, m);
    }

    // Reduce modulo X^p - X
    MultivariatePolynomial<piInt> remainder_poly{static_cast<std::size_t>(p), m};

    std::vector<std::size_t> decomposed_idx(m);
    for (const auto& coeff : poly.coeffs) {
        // Calculate index in new polynomial
        std::size_t new_idx = 0;
        std::size_t cur_factor = 1;
        for (std::size_t i = 0; i < decomposed_idx.size(); ++i) {
            if (decomposed_idx[i] <= static_cast<std::size_t>(p - 1)) {
                new_idx += decomposed_idx[i] * cur_factor;
            } else {
                new_idx += (decomposed_idx[i] % (p - 1)) * cur_factor;
            }
            cur_factor *= p;
        }

        remainder_poly.coeffs[new_idx] += coeff;

        // Increment index
        ++decomposed_idx[0];
        for (std::size_t j = 0; j < decomposed_idx.size() - 1; ++j) {
            if (decomposed_idx[j] == poly.d) {
                decomposed_idx[j] = 0;
                ++decomposed_idx[j + 1];
            } else {
                break;
            }
        }
    }
    return EvaluateAllImpl(remainder_poly.coeffs, m);
}

void MultidimensionalFFT::EvaluateUnivariate(std::span<const piInt> coeffs,
                                             std::span<piInt> out) const {
    std::vector<piInt> results(p - 1);
    if (p == 2) {
        results[0] = coeffs[0];
    } else {
        fft.transform(coeffs.data(), results.data());
    }

    int16_t cur = 1;
    for (std::size_t i = 0; i < results.size(); ++i) {
        out[cur] = results[i] + coeffs.back();
        cur = static_cast<int16_t>(static_cast<int32_t>(cur) * alpha % p);
    }
    out[0] = coeffs[0];
}

std::vector<piInt> MultidimensionalFFT::EvaluateAllImpl(std::span<const piInt> coeffs,
                                                        std::size_t stage) const {
    std::vector<piInt> out(coeffs.size());
    if (stage == 1) {
        EvaluateUnivariate(coeffs, out);
        return out;
    }
    std::vector<piInt> final_coeffs(coeffs.size());
    for (std::size_t i = 0; i < static_cast<std::size_t>(p); ++i) {
        auto evals =
            EvaluateAllImpl({coeffs.begin() + coeffs.size() / p * i, coeffs.size() / p}, stage - 1);
        for (std::size_t j = 0; j < evals.size(); ++j) {
            final_coeffs[j * p + i] = evals[j];
        }
    }

    for (std::size_t i = 0; i < coeffs.size() / p; ++i) {
        EvaluateUnivariate({final_coeffs.begin() + p * i, static_cast<std::size_t>(p)},
                           {out.begin() + p * i, static_cast<std::size_t>(p)});
    }
    return out;
}

// Use the segmented sieve of Eratosthenes
template <typename T>
static std::vector<T> GetPrimes(T N) {
    const auto delta = static_cast<std::size_t>(std::sqrt(N));

    std::vector<T> primes;
    std::vector<bool> flags(delta);
    // First segment
    for (std::size_t i = 2; i <= delta; ++i) {
        if (flags[i - 1]) {
            continue;
        }

        primes.emplace_back(static_cast<T>(i));
        for (std::size_t j = i * i; j <= delta; j += i) {
            flags[j - 1] = true;
        }
    }

    // Later segments
    for (std::size_t i = delta + 1; i <= static_cast<std::size_t>(N); i += delta) {
        flags.assign(delta, false);
        for (const auto prime : primes) {
            for (std::size_t j = static_cast<std::size_t>(((i - 1) / prime + 1) * prime);
                 j < i + delta && j <= static_cast<std::size_t>(N); j += prime) {

                flags[j - i] = true;
            }
            if (static_cast<std::size_t>(prime) * static_cast<std::size_t>(prime) >= i + delta) {
                break;
            }
        }
        for (std::size_t j = i; j < i + delta && j <= static_cast<std::size_t>(N); ++j) {
            if (!flags[j - i]) {
                primes.emplace_back(static_cast<T>(j));
            }
        }
    }
    return primes;
}

SimpleFastEvaluator::SimpleFastEvaluator(const std::filesystem::path& path_,
                                         const MultivariatePolynomial<pInt>& poly, int16_t p_)
    : p(p_), path(path_) {

    assert(p <= 16000);

    const auto logM = poly.m * std::log2(poly.d) + (poly.m * (poly.d - 1) + 1) * std::log2(p);
    primes = GetPrimes(static_cast<int16_t>(16 * logM + 1));

    {
        std::ofstream file{path / std::filesystem::u8path(PrimesFile), std::ios::binary};
        file.write(reinterpret_cast<const char*>(&p), sizeof(int16_t));

        const auto num_primes = primes.size();
        file.write(reinterpret_cast<const char*>(&num_primes), sizeof(std::size_t));
        file.write(reinterpret_cast<const char*>(primes.data()), primes.size() * sizeof(int16_t));
    }

    MultivariatePolynomial<piInt> fi(poly.d, poly.m);
    for (const auto pi : primes) {
        for (std::size_t j = 0; j < poly.coeffs.size(); ++j) {
            // Note that this reinterpretation is done in 0, ..., p - 1
            fi.coeffs[j] = piInt{static_cast<int16_t>((poly.coeffs[j].value + p) % p), pi};
        }

        MultidimensionalFFT fft{pi, poly.m};
        const auto& evaluation = fft.EvaluateAll(fi);

        std::ofstream file{path / std::filesystem::u8path(std::to_string(pi) + ".dat"),
                           std::ios::binary};
        for (const auto& it : evaluation) {
            file.write(reinterpret_cast<const char*>(&it.value), sizeof(int16_t));
        }
    }
}

SimpleFastEvaluator::SimpleFastEvaluator(const std::filesystem::path& path_) : path(path_) {
    std::ifstream file{path / std::filesystem::u8path(PrimesFile), std::ios::binary};
    file.read(reinterpret_cast<char*>(&p), sizeof(int16_t));

    std::size_t num_primes{};
    file.read(reinterpret_cast<char*>(&num_primes), sizeof(std::size_t));

    primes.resize(num_primes);
    file.read(reinterpret_cast<char*>(primes.data()), primes.size() * sizeof(int16_t));
}

SimpleFastEvaluator::~SimpleFastEvaluator() = default;

// Returns the result of the CRT modulo p
template <typename T>
static T ChineseRemainderTheorem(std::span<const T> results, std::span<const T> moduli, T p) {
    T cur_result = static_cast<T>(results[0] % p);
    T cur_modulus = static_cast<T>(moduli[0] % p);
    for (std::size_t i = 1; i < results.size(); ++i) {
        const auto result = static_cast<T>(results[i] % p);
        const auto modulus = static_cast<T>(moduli[i] % p);

        const auto& [_, m1, m2] = boost::integer::extended_euclidean(cur_modulus + p, modulus + p);

        // No need to implicitly promote for s16 but better safe than sorry
        using PromoteType =
            std::conditional_t<std::is_same_v<T, int16_t>, int32_t,
                               std::conditional_t<std::is_same_v<T, int32_t>, int64_t, T>>;

        cur_result = static_cast<T>((static_cast<PromoteType>(cur_result) * m2 % p * modulus % p +
                                     static_cast<PromoteType>(result) * m1 % p * cur_modulus % p) %
                                    p);
        cur_modulus = static_cast<T>(static_cast<PromoteType>(cur_modulus) * modulus % p);
    }
    return cur_result;
}

pInt SimpleFastEvaluator::Evaluate(std::span<const pInt> xs) const {
    std::vector<int16_t> crt_results;
    crt_results.reserve(primes.size());
    for (std::size_t i = 0; i < primes.size(); ++i) {
        const auto pi = primes[i];
        // Find the index in the table for that prime
        std::size_t idx = 0;
        std::size_t cur_factor = 1;
        for (auto it = xs.rbegin(); it != xs.rend(); ++it) {
            idx += ((it->value + p) % pi) * cur_factor;
            cur_factor *= pi;
        }

        std::ifstream file{path / std::filesystem::u8path(std::to_string(pi) + ".dat"),
                           std::ios::binary | std::ios::ate};
        const auto file_size = static_cast<std::size_t>(file.tellg());
        assert(file_size >= (idx + 1) * sizeof(int16_t));

        file.seekg(idx * sizeof(int16_t));

        int16_t value{};
        file.read(reinterpret_cast<char*>(&value), sizeof(int16_t));
        crt_results.emplace_back(value + pi);
    }
    return pInt{ChineseRemainderTheorem<int16_t>(std::span{crt_results}, std::span{primes}, p), p};
}

} // namespace Evaluation
