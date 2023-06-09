// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#pragma once

#include <filesystem>
#include <memory>
#include <span>
#include "core/common.h"
#include "core/kissfft.hpp"

namespace Evaluation {

// TODO: Move these tags into the classes, or make them templates or whatever, currently
// it's a little confusing

struct pIntTag {};
using pInt = ModularInt<int16_t, pIntTag>;
// p_i int
struct piIntTag {};
using piInt = ModularInt<int16_t, piIntTag>;

struct rIntTag {};
using rInt = ModularInt<BigInt, rIntTag>;

class FiniteFieldFFT;

class MultidimensionalFFT {
public:
    explicit MultidimensionalFFT(int16_t p, std::size_t m);
    ~MultidimensionalFFT();

    // Output go dimension 1 to m
    std::vector<piInt> EvaluateAll(const MultivariatePolynomial<piInt>& poly) const;

private:
    void EvaluateUnivariate(std::span<const piInt> coeffs, std::span<piInt> out) const;
    std::vector<piInt> EvaluateAllImpl(std::span<const piInt> coeffs, std::size_t stage) const;

    int16_t p{};
    std::size_t m{};
    int16_t alpha{}; // primitive root
    std::unique_ptr<FiniteFieldFFT> fft;
};

// Fast evaluator over Z_p
class SimpleFastEvaluator {
public:
    // Calculate from poly & p and save to folder.
    explicit SimpleFastEvaluator(std::filesystem::path path,
                                 const MultivariatePolynomial<pInt>& poly, int16_t p);
    // Load from folder. The files are machine-dependent (endianess)
    explicit SimpleFastEvaluator(std::filesystem::path path);
    ~SimpleFastEvaluator();

    pInt Evaluate(std::span<const pInt> xs) const;

private:
    std::filesystem::path path;
    int16_t p{};
    std::vector<int16_t> primes;

    static constexpr std::string_view PrimesFile = "primes.dat";
};

// Combines multiple direct evaluators to make an evaluator with log log r complexity
class ScalarFastEvaluator {
public:
    explicit ScalarFastEvaluator(std::filesystem::path path,
                                 const MultivariatePolynomial<rInt>& poly,
                                 std::shared_ptr<BigInt> r);
    explicit ScalarFastEvaluator(std::filesystem::path path);
    ~ScalarFastEvaluator();

    rInt Evaluate(std::span<const rInt> xs) const;

private:
    std::filesystem::path path;
    std::shared_ptr<BigInt> r;
    std::vector<int16_t> primes;
    std::vector<std::unique_ptr<SimpleFastEvaluator>> evaluators;

    static constexpr std::string_view PrimesFile = "primes.dat";
};

// E(Y) = Y^e + 1
class UnivariateFastEvaluator {
public:
    using Element = UnivariateRingPolynomial<qInt>;

    explicit UnivariateFastEvaluator(std::filesystem::path path,
                                     const MultivariatePolynomial<Element>& poly,
                                     std::shared_ptr<BigInt> q, std::size_t e);
    explicit UnivariateFastEvaluator(std::filesystem::path path);
    ~UnivariateFastEvaluator();

    Element Evaluate(std::span<const Element> xs) const;

private:
    std::filesystem::path path;
    std::shared_ptr<BigInt> q;
    std::size_t e{};
    BigInt M;
    std::size_t D{};
    std::shared_ptr<BigInt> r;
    std::unique_ptr<ScalarFastEvaluator> evaluator;

    static constexpr std::string_view MetadataFile = "univariate.dat";
};

// E(Y) = Y^e1 + 1, E(Z) = Z^e2 + 1
class BivariateFastEvaluator {
public:
    using Element = BivariateRingPolynomial<qInt>;

    explicit BivariateFastEvaluator(std::filesystem::path path,
                                    const MultivariatePolynomial<Element>& poly,
                                    std::shared_ptr<BigInt> q, std::size_t e1, std::size_t e2);
    explicit BivariateFastEvaluator(std::filesystem::path path);
    ~BivariateFastEvaluator();

    Element Evaluate(std::span<const Element> xs) const;

private:
    std::filesystem::path path;
    std::shared_ptr<BigInt> q;
    std::size_t e1{};
    std::size_t e2{};
    BigInt M;
    std::size_t D1{};
    std::size_t D2{};
    std::shared_ptr<BigInt> r;
    std::unique_ptr<UnivariateFastEvaluator> evaluator;

    static constexpr std::string_view MetadataFile = "bivariate.dat";
};

} // namespace Evaluation
