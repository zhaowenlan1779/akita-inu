// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#include <catch2/catch_test_macros.hpp>
#include "core/big_int.h"

TEST_CASE("big int arithmetic", "[base]") {
    SECTION("small") {
        BigInt a = 5;
        BigInt b = 7;
        REQUIRE(a != b);
        REQUIRE(a + b == 12);
        REQUIRE(a - b == -2);
        REQUIRE(a * b == 35);
    }

    SECTION("large") {
        BigInt a{"190090923423049230423423"};
        BigInt b{"98934892384923476856"};
        REQUIRE(a != b);
        REQUIRE(a + b == BigInt{"190189858315434153900279"});
        REQUIRE(b - a == BigInt{"-189991988530664306946567"});
        REQUIRE(a * b == BigInt{"18806625052210105085392227519155074220798088"});
    }
}
