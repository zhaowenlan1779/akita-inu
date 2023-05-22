// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include <span>
#include "core/common.h"
#include "core/kissfft.hpp"

namespace Evaluation {

class MultidimensionalFFT {
public:
    explicit MultidimensionalFFT(int64_t d, std::size_t m);
    ~MultidimensionalFFT();

    // Output go dimension 1 to m
    std::vector<dInt> EvaluateAll(const MultiPoly& poly) const;

private:
    void EvaluateUnivariate(std::span<const dInt> coeffs, std::span<dInt> out) const;
    std::vector<dInt> EvaluateAllImpl(std::span<const dInt> coeffs, std::size_t m) const;

    int64_t d{};
    std::size_t m{};
    int64_t alpha{}; // primitive root
    kissfft<dInt, dInt> fft;
};

} // namespace Evaluation
