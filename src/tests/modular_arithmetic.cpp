// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <catch2/catch_test_macros.hpp>
#include "core/modular_arithmetic.hpp"

TEST_CASE("modular arithmetic", "[base]") {
    SECTION("odd modulus") {
        ModularInt<int> a{88, 97};
        ModularInt<int> b{17, 97};

        REQUIRE(a != b);
        REQUIRE(a + b == ModularInt<int>{8, 97});
        REQUIRE(b - a == ModularInt<int>{26, 97});
        REQUIRE(3 - b == ModularInt<int>{-14, 97});
        REQUIRE(a.GetRepresentative() == 39);
        REQUIRE(a.Reinterpret(17) == ModularInt<int>{5, 17});

        ModularInt<int> c{3, 97};
        REQUIRE(c.Reinterpret(17) == ModularInt<int>{5, 17});
    }

    SECTION("even modulus reinterpret") {
        ModularInt<int> a{88, 96};
        REQUIRE(a.GetRepresentative() == 40);
        REQUIRE(a.Reinterpret(17) == ModularInt<int>{6, 17});

        ModularInt<int> c{4, 96};
        REQUIRE(c.Reinterpret(17) == ModularInt<int>{7, 17});
    }
}
