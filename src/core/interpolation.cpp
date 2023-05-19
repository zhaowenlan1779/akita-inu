// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <algorithm>
#include "core/interpolation.h"

UnivariateInterpolator::UnivariateInterpolator(int64_t d_) : d(d_) {
    Build(1, 0, d - 1);
}

UnivariateInterpolator::~UnivariateInterpolator() = default;

static std::pair<UnivariateInterpolator::dInt, std::size_t> FindHighestTerm(
    const UnivariateInterpolator::Poly& poly) {
    const auto it = std::find_if(poly.coeffs.rbegin(), poly.coeffs.rend(),
                                 [](const auto& x) { return static_cast<bool>(x); });
    return {*it, poly.coeffs.rend() - it - 1};
}

UnivariateInterpolator::Poly UnivariateInterpolator::Remainder(Poly a, Poly b) const {
    const auto& [coeff_b, deg_b] = FindHighestTerm(b);
    while (true) {
        const auto& [coeff_a, deg_a] = FindHighestTerm(a);
        if (deg_a < deg_b) {
            return a;
        }

        Poly monic{a.n};
        monic.coeffs[deg_a - deg_b] = dInt{1, d};
        a -= coeff_a / coeff_b * b * monic;
    }
}

void UnivariateInterpolator::Build(std::size_t k, std::size_t l, std::size_t r) {
    if (l == r) {
        Poly m{static_cast<std::size_t>(d + 1)};
        m.coeffs[0] = -dInt{static_cast<int64_t>(l), d};
        m.coeffs[1] = dInt{1, d};

        if (M.size() <= k) {
            M.resize(k + 1);
        }
        M[k] = std::move(m);
        return;
    }

    const auto mid = (l + r) / 2;
    Build(2 * k, l, mid);
    Build(2 * k + 1, mid + 1, r);
    M[k] = M[2 * k] * M[2 * k + 1];
}

std::vector<UnivariateInterpolator::dInt> UnivariateInterpolator::Evaluate(Poly poly, std::size_t k,
                                                                           std::size_t l,
                                                                           std::size_t r) const {

    if (l == r) {
        return {std::move(poly.coeffs[0])};
    }

    const auto mid = (l + r) / 2;
    auto r_0 = Remainder(poly, M[2 * k]);
    auto result = Evaluate(std::move(r_0), 2 * k, l, mid);

    auto r_1 = Remainder(std::move(poly), M[2 * k + 1]);
    auto right = Evaluate(std::move(r_1), 2 * k + 1, mid + 1, r);

    result.resize(result.size() + right.size());
    std::move(right.begin(), right.end(), result.end() - right.size());

    return result;
}

UnivariateInterpolator::Poly UnivariateInterpolator::CombineModuli(std::span<dInt> c, std::size_t k,
                                                                   std::size_t l,
                                                                   std::size_t r) const {
    if (l == r) {
        Poly poly{static_cast<std::size_t>(d + 1)};
        poly.coeffs[0] = std::move(c[0]);
        return poly;
    }

    const auto mid = (l + r) / 2;
    auto r_0 = CombineModuli({c.begin(), (c.size() + 1) / 2}, 2 * k, l, mid);
    auto r_1 = CombineModuli({c.begin() + (c.size() + 1) / 2, c.end()}, 2 * k + 1, mid + 1, r);
    return M[2 * k + 1] * r_0 + M[2 * k] * r_1;
}

UnivariateInterpolator::Poly UnivariateInterpolator::Interpolate(std::span<dInt> ys) const {
    auto m = M[1];

    // Take derivative of m
    for (std::size_t i = 0; i < m.coeffs.size() - 1; ++i) {
        m.coeffs[i] = m.coeffs[i + 1] * dInt{static_cast<int64_t>(i + 1), d};
    }
    m.coeffs.back() = dInt{0, d};

    // Evaluate m'
    const auto& s_inv = Evaluate(std::move(m), 1, 0, d - 1);

    std::vector<dInt> c(s_inv.size());
    for (std::size_t i = 0; i < s_inv.size(); ++i) {
        c[i] = ys[i] / s_inv[i];
    }
    return CombineModuli(c, 1, 0, d - 1);
}
