// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include <span>
#include "core/modular_arithmetic.hpp"
#include "core/ring_polynomial.hpp"

// Algorithm from Joachim von zur Gathen and Jurgen Gerhard. Modern computer algebra. Cambridge
// University Press, Cambridge, England, 3rd edition, 2013. Chapter 10.
// This version interpolates on fixed points 0, 1, ..., d - 1
class UnivariateInterpolator {
public:
    struct dIntTag {};
    using dInt = ModularInt<int64_t, dIntTag, true>;

    using Poly = UnivariatePolynomial<dInt>;

    explicit UnivariateInterpolator(int64_t d);
    ~UnivariateInterpolator();

    Poly Interpolate(std::span<dInt> ys) const;

private:
    Poly Remainder(Poly a, Poly b) const;

    void Build(std::size_t k, std::size_t l, std::size_t r);
    // Recursive algorithm from Algorithm 10.5.
    std::vector<dInt> Evaluate(Poly poly, std::size_t k, std::size_t l, std::size_t r) const;
    Poly CombineModuli(std::span<dInt> c, std::size_t k, std::size_t l, std::size_t r) const;

    // modulus
    int64_t d{};
    // subproduct tree. This is just like a segment tree.
    std::vector<Poly> M;
};
