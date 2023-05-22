// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <cmath>
#include <iostream>
#include <boost/integer/common_factor.hpp>
#include "core/ASHE.h"
#include "core/kissfft.hpp"

ASHE::ASHE(std::size_t lambda, int64_t d_, std::size_t D_, std::size_t N) : D(D_), d(d_) {
    // Setup parameters
    const auto logd = static_cast<std::size_t>(std::ceil(std::log2(d)));
    const auto t = D * logd + std::log2(N) + logd + 1;

    // Raise to next 2's power. This does not make it very much bigger.
    auto target_n = static_cast<std::size_t>(lambda + t);
    n = 1;
    while (target_n >>= 1) {
        n <<= 1;
    }
    n <<= 1;

    // Select q as smallest integer larger than (2n^(3/2) log n)^t coprime to d
    BigInt tentative_q = boost::multiprecision::pow(
        BigInt{std::ceil(2 * n * std::sqrt(n) * std::log2(n))}, static_cast<unsigned int>(t) + 1);
    while (boost::math::gcd(tentative_q, d) != 1) {
        ++tentative_q;
    }
    q = std::make_shared<BigInt>(std::move(tentative_q));
}

ASHE::~ASHE() = default;

ASHE::QElement ASHE::Gen() const {
    return SampleUniform();
}

ASHE::RElement ASHE::Enc(const QElement& key, const dInt& plaintext) const {
    const auto a = SampleUniform();
    const auto b = a * key + qInt{d, q} * SampleError() + plaintext.Reinterpret<qInt>(q);

    // Construct ct = -aY + b
    RElement ct{n, D};
    for (std::size_t i = 0; i < n; ++i) {
        ct.Coeff(i, 0) = b.coeffs[i];
        ct.Coeff(i, 1) = -a.coeffs[i];
    }
    return ct;
}

dInt ASHE::Dec(const QElement& key, const RElement& ct) const {
    // Evaluate ct on key
    QElement result{n};

    QElement cur{n};
    cur.coeffs[0] = qInt{1, q};

    for (std::size_t j = 0; j < D; ++j) {
        // Reconstruct degree j as QElement
        QElement k{n};
        for (std::size_t i = 0; i < n; ++i) {
            k.coeffs[i] = ct.Coeff(i, j);
        }

        result += k * cur;
        if (j != D - 1) {
            cur *= key;
        }
    }

    return result.coeffs[0].Reinterpret<dInt>(d);
}

ASHE::RElement ASHE::Lift(const dInt& val) const {
    RElement ct{n, D};
    ct.Coeff(0, 0) = val.Reinterpret<qInt>(q);
    return ct;
}

ASHE::QElement ASHE::SampleUniform() const {
    boost::random::uniform_int_distribution<BigInt> dist{0, *q - 1};

    QElement out{n};
    for (std::size_t i = 0; i < n; ++i) {
        out.coeffs[i] = qInt{dist(gen), q};
    }
    return out;
}

constexpr double PI = 3.141592653589793238462643383279502884L;

bool ASHE::SampleErrorTrial(QElement& out) const {
    // Reference:
    // Chris Peikert, Oded Regev, and Noah Stephens-Davidowitz. Pseudorandomness
    // of ring-LWE for any ring and modulus. In Hamed Hatami, Pierre McKenzie, and
    // Valerie King, editors, 49th Annual ACM Symposium on Theory of Computing, pages
    // 461-473. ACM Press, June 2017.

    // Selected values: alpha = 2 / q
    const auto f2 = std::log2(n);

    // The standard deviation is sqrt(n log n), we add another log factor to make it
    // statistically close as n approaches infinity
    const auto beta = std::sqrt(n) * std::log2(n);

    // First sample a distribution from the family.
    // Since n is even, there are no real embeddings, and n/2 pairs of complex embeddings
    std::vector<double> r(n / 2);
    for (auto& it : r) {
        boost::random::normal_distribution dist{0.0, 1 / 2 / std::sqrt(PI)};

        const auto x = dist(gen), y = dist(gen);
        // For numerical stability alpha2 is omitted here
        it = std::sqrt((x * x + y * y + f2) / 2);
    }

    // Sample an embedding from the distribution.
    // The values are 0 for even indices, such that inverse DFT can recover the coefficients
    // from the embedding
    std::vector<std::complex<double>> embedding(2 * n);
    for (std::size_t i = 0; i < r.size(); ++i) {
        boost::random::normal_distribution dist{0.0, r[i] / 2 / std::sqrt(PI)};

        using namespace std::complex_literals;
        embedding[2 * i + 1] = dist(gen) + dist(gen) * 1i;
        embedding[2 * (n - 1 - i) + 1] = std::conj(embedding[2 * i + 1]);
    }

    // Inverse DFT on the embeddings to obtain coefficients * n / alpha
    kissfft<double> fft{embedding.size(), true};

    std::vector<std::complex<double>> coeffs(2 * n);
    fft.transform(embedding.data(), coeffs.data());
    embedding.clear();

    for (std::size_t i = 0; i < n; ++i) {
        // alpha = 2 / q => q * r * alpha / n = 2 * r / n
        // We multiply this result by n because we are working in R and not R^\vee.
        // Ref: A Toolkit for Ring-LWE Cryptography, Page 6; Ring-LWE in Polynomial Rings, Page 39
        const auto value = 2 * coeffs[i].real();

        if (std::abs(std::round(value)) > beta) {
            return false;
        }
        out.coeffs[i] = qInt{static_cast<int64_t>(std::round(value)), q};
    }
    return true;
}

ASHE::QElement ASHE::SampleError() const {
    QElement result{n};
    std::size_t i = 0;
    while (!SampleErrorTrial(result)) {
        std::cout << "Warning: Sample error failed. Retrying (Attempt " << ++i << ")." << std::endl;
    }
    return result;
}
