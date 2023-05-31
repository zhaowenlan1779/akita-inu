// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#include <cmath>
#include <ranges>
#include "core/DEPIR.h"
#include "core/interpolation.h"

namespace DEPIR {

Server::Server(const std::span<dInt>& raw_db, int64_t d_) : d(d_) {
    m = static_cast<std::size_t>(std::ceil(std::log2(raw_db.size()) / std::log2(d)));

    Interpolation::MultivariateInterpolator interp{d, m};
    db = interp.Interpolate(raw_db);
}

Server::~Server() = default;

std::pair<int64_t, std::size_t> Server::GetServerParams() const {
    return {d, m};
}

void Server::SetClientParams(std::size_t n_, std::size_t D_, BigInt q_) {
    n = n_;
    D = D_;

    q = std::make_shared<BigInt>();
    *q = std::move(q_);
}

ASHE::RElement Server::Resp(std::span<const ASHE::RElement> ct_orig) const {
    // Evaluate db slowly
    // While there is an implementation of the fast evaluator, it takes up so much storage and
    // computation that it's practically unusable

    // The order of dimensions is reversed in interp.Interpolate. Hence we need to reverse here.
    const auto& ct = ct_orig | std::views::reverse;

    // Current powers of the variables
    std::vector<ASHE::RElement> cur_ct(ct.size());
    for (auto& it : cur_ct) {
        it = ASHE::RElement{n, D};
        it.coeffs[0] = qInt{1, q};
    }

    std::vector<std::size_t> decomposed_idx(m);
    ASHE::RElement result{n, D};
    for (const auto& coeff : db.coeffs) {
        ASHE::RElement cur_term{n, D};
        cur_term.coeffs[0] = coeff.Reinterpret<qInt>(q);
        for (const auto& it : cur_ct) {
            cur_term *= it;
        }
        result += cur_term;

        // Increment index and update the current powers
        for (std::size_t j = 0; j < decomposed_idx.size(); ++j) {
            ++decomposed_idx[j];
            if (decomposed_idx[j] == static_cast<std::size_t>(d)) {
                decomposed_idx[j] = 0;
                cur_ct[j] = ASHE::RElement{n, D};
                cur_ct[j].coeffs[0] = qInt{1, q};
            } else {
                cur_ct[j] *= ct[j];
                break;
            }
        }
    }
    return result;
}

Client::Client(std::size_t lambda, int64_t d_, std::size_t m_)
    : d(d_), m(m_), ashe(lambda, d, d * m, static_cast<std::size_t>(std::pow(d, m))) {}

Client::~Client() = default;

std::tuple<std::size_t, std::size_t, BigInt> Client::GetClientParams() const {
    return {ashe.n, ashe.D, *ashe.q};
}

std::pair<std::vector<ASHE::RElement>, ASHE::QElement> Client::Query(std::size_t idx) const {
    const auto s = ashe.Gen();
    std::vector<ASHE::RElement> ct(m);
    for (std::size_t i = 0; i < m; ++i) {
        ct[i] = ashe.Enc(s, dInt{static_cast<int64_t>(idx) % d, d});
        idx /= d;
    }
    return {ct, s};
}

dInt Client::Dec(const ASHE::QElement& s, const ASHE::RElement& ans) const {
    return ashe.Dec(s, ans);
}

} // namespace DEPIR
