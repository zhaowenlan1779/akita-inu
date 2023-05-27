// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#include <catch2/catch_test_macros.hpp>
#include "core/modular_arithmetic.hpp"

TEST_CASE("modular arithmetic", "[base]") {
    SECTION("odd modulus") {
        ModularInt<int> a{88, 97};
        ModularInt<int> b{17, 97};

        REQUIRE(a != b);
        REQUIRE(a + b == ModularInt<int>{8, 97});
        REQUIRE(b - a == ModularInt<int>{26, 97});
        REQUIRE(a.GetRepresentative() == -9);

        REQUIRE(a.Reinterpret(17) == ModularInt<int>{8, 17});

        ModularInt<int> c{3, 97};
        REQUIRE(c.Reinterpret(17) == ModularInt<int>{3, 17});
    }

    SECTION("even modulus reinterpret") {
        ModularInt<int> a{88, 96};
        REQUIRE(a.GetRepresentative() == -8);
        REQUIRE(a.Reinterpret(17) == ModularInt<int>{9, 17});

        ModularInt<int> c{4, 96};
        REQUIRE(c.Reinterpret(17) == ModularInt<int>{4, 17});
    }

    SECTION("default zeros") {
        ModularInt<int> def;
        ModularInt<int> a{88, 96};

        REQUIRE(a + def == a);
        REQUIRE(a - def == a);
        REQUIRE(a * def == ModularInt<int>{0, 96});
        REQUIRE(def + a == a);
        REQUIRE(def - a == -a);
        REQUIRE(def * a == ModularInt<int>{0, 96});
    }
}
