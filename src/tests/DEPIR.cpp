// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#include <iostream>
#include <catch2/catch_test_macros.hpp>
#include "core/DEPIR.h"

TEST_CASE("DEPIR", "[DEPIR]") {
    static constexpr int64_t d = 3;

    std::vector<dInt> raw_db(3 * 3 * 3);
    for (auto& it : raw_db) {
        it = dInt{std::rand() % d, d};
    }

    // Set up parameters
    DEPIR::Server server{raw_db, d};
    const auto& [_, m] = server.GetServerParams();
    DEPIR::Client client{0, d, m};
    const auto& [n, D, q] = client.GetClientParams();
    server.SetClientParams(n, D, q);

    for (std::size_t idx = 0; idx < raw_db.size(); ++idx) {
        const auto& [query, key] = client.Query(idx);
        const auto& resp = server.Resp(query);
        REQUIRE(raw_db[idx] == client.Dec(key, resp));
    }
}
