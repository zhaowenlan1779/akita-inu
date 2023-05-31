// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#include <span>
#include <tuple>
#include "core/ASHE.h"
#include "core/common.h"

namespace DEPIR {

class Server {
public:
    // The size of the db must be a power of d
    explicit Server(const std::span<dInt>& raw_db, int64_t d);
    ~Server();

    // (d, m)
    std::pair<int64_t, std::size_t> GetServerParams() const;
    void SetClientParams(std::size_t n, std::size_t D, BigInt q);

    ASHE::RElement Resp(std::span<const ASHE::RElement> ct) const;

private:
    int64_t d{};
    std::size_t m{};
    std::size_t n{};
    std::size_t D{};
    std::shared_ptr<BigInt> q{};
    MultiPoly db;
};

class Client {
public:
    explicit Client(std::size_t lambda, int64_t d, std::size_t m);
    ~Client();

    // (n, D, q)
    std::tuple<std::size_t, std::size_t, BigInt> GetClientParams() const;

    std::pair<std::vector<ASHE::RElement>, ASHE::QElement> Query(std::size_t idx) const;
    dInt Dec(const ASHE::QElement& s, const ASHE::RElement& ans) const;

private:
    int64_t d{};
    std::size_t m{};
    ASHE ashe;
};

} // namespace DEPIR
