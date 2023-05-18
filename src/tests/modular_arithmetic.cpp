// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <catch2/catch_test_macros.hpp>
#include "core/modular_arithmetic.hpp"

TEST_CASE("modular arithmetic", "[base]") {
    SECTION("odd modulus") {
        const auto mod97 = std::make_shared<int>(97);
        ModularInt<int> a{88, mod97};
        ModularInt<int> b{17, mod97};

        REQUIRE(a != b);
        REQUIRE(a + b == ModularInt<int>{8, mod97});
        REQUIRE(b - a == ModularInt<int>{26, mod97});
        REQUIRE(a.GetRepresentative() == -9);

        const auto mod17 = std::make_shared<int>(17);
        REQUIRE(a.Reinterpret(mod17) == ModularInt<int>{8, mod17});

        ModularInt<int> c{3, mod97};
        REQUIRE(c.Reinterpret(mod17) == ModularInt<int>{3, mod17});
    }

    SECTION("even modulus reinterpret") {
        const auto mod96 = std::make_shared<int>(96);

        ModularInt<int> a{88, mod96};
        REQUIRE(a.GetRepresentative() == -8);

        const auto mod17 = std::make_shared<int>(17);
        REQUIRE(a.Reinterpret(mod17) == ModularInt<int>{9, mod17});

        ModularInt<int> c{4, mod96};
        REQUIRE(c.Reinterpret(mod17) == ModularInt<int>{4, mod17});
    }

    SECTION("default zeros") {
        ModularInt<int> def;

        const auto mod96 = std::make_shared<int>(96);
        ModularInt<int> a{88, mod96};

        REQUIRE(a != def);
        REQUIRE(a + def == a);
        REQUIRE(a - def == a);
        REQUIRE(a * def == ModularInt<int>{0, mod96});
        REQUIRE(def + a == a);
        REQUIRE(def - a == -a);
        REQUIRE(def * a == ModularInt<int>{0, mod96});
    }
}
