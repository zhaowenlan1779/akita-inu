// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#include <algorithm>
#include <fstream>
#include <type_traits>
#include <BS_thread_pool_light.hpp>
#include <boost/integer/extended_euclidean.hpp>
#include "core/fast_evaluation.h"
#include "core/kissfft.hpp"

namespace Evaluation {

// Sequences in F_p of length p-1, forward DFT
class FiniteFieldFFT {
public:
    virtual ~FiniteFieldFFT() = default;
    virtual void Transform(std::span<const piInt> src, std::span<piInt> dst) const = 0;
};

namespace {

class FiniteFieldFFTDirect : public FiniteFieldFFT {
public:
    explicit FiniteFieldFFTDirect(int16_t p, int16_t alpha);
    ~FiniteFieldFFTDirect() override = default;
    void Transform(std::span<const piInt> src, std::span<piInt> dst) const override;

private:
    int16_t p{};
    kissfft<double, false> fft;
    mutable std::vector<std::complex<double>> src_scratch;
    mutable std::vector<std::complex<double>> dst_scratch;
};

FiniteFieldFFTDirect::FiniteFieldFFTDirect(int16_t p_, int16_t alpha)
    : p(p_), fft{static_cast<std::size_t>(p - 1), false}, src_scratch(p - 1), dst_scratch(p - 1) {
    if (p == 2) {
        return;
    }

    // Initialize twiddles
    fft._twiddles.resize(p - 1);
    piInt cur{1, p};
    for (int16_t i = 0; i < p - 1; ++i) {
        fft._twiddles[i] = cur.value;
        cur *= alpha;
    }
}

void FiniteFieldFFTDirect::Transform(std::span<const piInt> src, std::span<piInt> dst) const {
    if (p == 2) {
        dst[0] = src[0];
        return;
    }

    for (std::size_t i = 0; i < src.size(); ++i) {
        src_scratch[i] = src[i].value;
    }
    fft.transform(src_scratch.data(), dst_scratch.data());
    for (std::size_t i = 0; i < dst.size(); ++i) {
        dst[i] = {static_cast<int16_t>(static_cast<int64_t>(std::round(dst_scratch[i].real())) % p),
                  p};
    }
}

// Uses Bluestein's algorithm, and some techniques described in The Fast Fourier Transform in a
// Finite Field, by J. M. Pollard.
class FiniteFieldFFTBluestein : public FiniteFieldFFT {
public:
    explicit FiniteFieldFFTBluestein(int16_t p, int16_t alpha);
    ~FiniteFieldFFTBluestein() override = default;
    void Transform(std::span<const piInt> src, std::span<piInt> dst) const override;

private:
    void Convolve(std::span<const double> a, std::span<piInt> out) const;

    int16_t p{};
    int16_t p1{};
    kissfft<double> fft;
    kissfft<double> inv_fft;
    std::vector<piInt> twiddles;
    std::vector<std::complex<double>> s_dft;

    // scratch buffers, to avoid dynamic allocation
    mutable std::vector<double> even;
    mutable std::vector<double> odd;
    mutable std::vector<std::complex<double>> scratch1;
    mutable std::vector<std::complex<double>> scratch2;
    mutable std::vector<piInt> out_even;
    mutable std::vector<piInt> out_odd;
};

// These are the primes p where p-1 is highly composite (with only prime factors 2, 3, 5)
constexpr std::array<int16_t, 65> GoodPrimesList{{
    2,    3,    5,    7,    11,   13,   17,   19,    31,    37,    41,    61,    73,
    97,   101,  109,  151,  163,  181,  193,  241,   251,   257,   271,   401,   433,
    487,  541,  577,  601,  641,  751,  769,  811,   1153,  1201,  1297,  1459,  1601,
    1621, 1801, 2161, 2251, 2593, 2917, 3001, 3457,  3889,  4001,  4051,  4801,  4861,
    6481, 7681, 8101, 8641, 9001, 9601, 9721, 10369, 11251, 12289, 14401, 15361, 16001,
}};

int16_t FindNextGoodPrime(int16_t n) {
    assert(n <= 16000);
    return *std::lower_bound(GoodPrimesList.begin(), GoodPrimesList.end(), n);
}

FiniteFieldFFTBluestein::FiniteFieldFFTBluestein(int16_t p_, int16_t alpha)
    : p(p_), p1(FindNextGoodPrime(p)), fft((p1 - 1) / 2, false), inv_fft(p1 - 1, true),
      s_dft((p1 - 1) / 2), even(p1 - 1), odd(p1 - 1), scratch1(p1 - 1), scratch2(p1 - 1),
      out_even(p1 - 1), out_odd(p1 - 1) {

    if (p == 2) {
        return;
    }

    twiddles.resize(p - 1);
    piInt cur{1, p};
    for (int16_t i = 0; i < p - 1; ++i) {
        twiddles[i] = cur;
        cur *= alpha;
    }

    // Generate s (used for convolution)
    std::vector<double> s(p1 - 1);
    for (std::size_t i = 0; i < static_cast<std::size_t>((p - 1) / 2); ++i) {
        // s^{-i^2}
        const auto idx = (static_cast<int64_t>(twiddles.size()) - static_cast<int64_t>(i * i) +
                          static_cast<int64_t>(p - 1) * (p - 1)) %
                         (p - 1);
        s[i] = twiddles[idx].GetRepresentative();
        if (i > 0) {
            s[(p1 - 1) - i] = s[i];
        }
    }
    fft.transform_real(s.data(), s_dft.data());
}

void FiniteFieldFFTBluestein::Convolve(std::span<const double> a, std::span<piInt> out) const {
    fft.transform_real(a.data(), scratch1.data());

    // See kissfft doc for interpretation of the values
    const auto a0 = scratch1[0];
    scratch1[0] = a0.real() * s_dft[0].real();
    scratch1[(p1 - 1) / 2] = a0.imag() * s_dft[0].imag();

    for (std::size_t i = 1; i < scratch1.size() / 2; ++i) {
        scratch1[i] *= s_dft[i];
    }
    for (std::size_t i = 1; i < scratch1.size() / 2; ++i) {
        scratch1[scratch1.size() - i] = std::conj(scratch1[i]);
    }

    inv_fft.transform(scratch1.data(), scratch2.data());

    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = {static_cast<int16_t>(
                      static_cast<int64_t>(std::round(scratch2[i].real() / scratch1.size())) % p),
                  p};
        out[i] *= twiddles[(i * i) % (p - 1)];
    }
}

void FiniteFieldFFTBluestein::Transform(std::span<const piInt> src, std::span<piInt> dst) const {
    if (p == 2) {
        dst[0] = src[0];
        return;
    }

    // Pick out even and odd terms
    for (std::size_t i = 0; i < src.size() / 2; ++i) {
        even[i] = (src[2 * i] * twiddles[(i * i) % (p - 1)]).GetRepresentative();
        odd[i] = (src[2 * i + 1] * twiddles[(i * i) % (p - 1)]).GetRepresentative();
    }

    Convolve(even, out_even);
    Convolve(odd, out_odd);
    for (std::size_t i = 0; i < src.size() / 2; ++i) {
        dst[i] = out_even[i] + twiddles[i] * out_odd[i];
        dst[i + src.size() / 2] = out_even[i] - twiddles[i] * out_odd[i];
    }
}

double EstimateCost(int16_t n) {
    static constexpr double Penalty = 1.1; // penalty for non-hardcoded larger factors

    const auto original_n = n;

    double result = 0;
    while ((n & 1) == 0) {
        result += 2;
        n >>= 1;
    }
    for (std::size_t x = 3; x * x <= static_cast<std::size_t>(n); x += 2) {
        while (n % x == 0) {
            // penalize larger prime factors
            result += (x <= 5) ? double(x) : Penalty * double(x);
            n /= x;
        }
    }
    if (n > 1)
        result += (n <= 5) ? double(n) : Penalty * double(n);
    return result * original_n;
}

std::unique_ptr<FiniteFieldFFT> CreateFiniteFieldFFT(int16_t p, int16_t alpha) {
    // Estimate cost
    // This factor of 2 was verified by experiment
    const auto cost_direct = EstimateCost(p - 1);
    const auto cost_bluestein = 2 * EstimateCost(FindNextGoodPrime(p) - 1);

    if (cost_direct <= cost_bluestein) {
        return std::make_unique<FiniteFieldFFTDirect>(p, alpha);
    } else {
        return std::make_unique<FiniteFieldFFTBluestein>(p, alpha);
    }
}

std::vector<int16_t> GetDistinctPrimeFactors(int16_t s) {
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
            if (static_cast<int32_t>(p) * p > s) {
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

int16_t QuickPow(int16_t base, int16_t power, int16_t modulus) {
    int32_t result = 1;
    int32_t cur_base = base;
    while (power > 0) {
        if (power & 1) {
            result *= cur_base;
            result %= modulus;
        }
        power >>= 1;
        cur_base *= cur_base;
        cur_base %= modulus;
    }
    return static_cast<int16_t>(result);
}

int16_t GetPrimitiveRoot(int16_t p) {
    if (p == 2) {
        return 2; // irrelevant
    }

    const auto s = p - 1;
    const auto& factors = GetDistinctPrimeFactors(s);

    // Find the smallest primitive root
    for (int16_t t = 2; t < p; ++t) {
        if (std::all_of(factors.begin(), factors.end(),
                        [p, t, s](int16_t q) { return QuickPow(t, s / q, p) != 1; })) {
            return t;
        }
    }
    throw std::runtime_error("Failed to find primitive root");
}

} // namespace

MultidimensionalFFT::MultidimensionalFFT(int16_t p_, std::size_t m_)
    : p(p_), m(m_), alpha(GetPrimitiveRoot(p)) {

    fft = CreateFiniteFieldFFT(p, alpha);
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
    fft->Transform({coeffs.begin(), coeffs.begin() + p - 1}, results);

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

    BS::thread_pool_light pool;
    for (const auto pi : primes) {
        const auto task = [this, pi, &poly] {
            MultivariatePolynomial<piInt> fi(poly.d, poly.m);
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
        };
        if (pi <= 100) {
            task();
        } else {
            pool.push_task(task);
        }
    }
    pool.wait_for_tasks();
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

        // No need to explicitly promote for s16 but better safe than sorry
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
