// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include <boost/random.hpp>
#include "core/big_int.h"
#include "core/common.h"
#include "core/modular_arithmetic.hpp"
#include "core/ring_polynomial.hpp"

class ASHE {
public:
    struct qIntTag {};
    using qInt = ModularInt<BigInt, qIntTag>;

    using QElement = UnivariateRingPolynomial<qInt>;
    using RElement = BivariateRingPolynomial<qInt>;

    explicit ASHE(std::size_t lambda, int64_t d, std::size_t D, std::size_t N);
    ~ASHE();

    QElement Gen() const;
    RElement Enc(const QElement& key, const dInt& plaintext) const;
    dInt Dec(const QElement& key, const RElement& ct) const;
    RElement Lift(const dInt& val) const;

    // Public parameters
    std::size_t n{}, D{};
    int64_t d{};
    std::shared_ptr<BigInt> q;

private:
    QElement SampleUniform() const;
    bool SampleErrorTrial(QElement& out) const;
    QElement SampleError() const;

    mutable boost::random::mt19937 gen;
};
