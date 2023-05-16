// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#include <catch2/catch_test_macros.hpp>
#include "core/ring_polynomial.hpp"

TEST_CASE("univariate ring polynomial", "[base]") {
    UnivariateRingPolynomial<int> a{3};
    a.coeffs = {7, 9, 12}; // 12z^2 + 9z + 7

    UnivariateRingPolynomial<int> b{3};
    b.coeffs = {-5, 3, -15}; // -15z^2 + 3z - 5

    REQUIRE(a != b);

    {
        UnivariateRingPolynomial<int> c{3};
        c.coeffs = {2, 12, -3};
        REQUIRE(a + b == c);
    }

    {
        UnivariateRingPolynomial<int> c{3};
        c.coeffs = {12, 6, 27};
        REQUIRE(a - b == c);
    }

    {
        UnivariateRingPolynomial<int> c{3};
        c.coeffs = {64, 156, -138};
        REQUIRE(a * b == c);
    }
}

TEST_CASE("bivariate ring polynomial", "[base]") {
    BivariateRingPolynomial<int> a{3, 2};
    a.coeffs = {7, 9, 12, 3, 4, 5}; // 7 + 9y + 12z + 3yz + 4z^2 + 5yz^2

    BivariateRingPolynomial<int> b{3, 2};
    b.coeffs = {-5, 3, -15, 1, 2, 3}; // -5 + 3y -15z + yz + 2z^2 + 3yz^2

    REQUIRE(a != b);

    {
        BivariateRingPolynomial<int> c{3, 2};
        c.coeffs = {2, 12, -3, 4, 6, 8};
        REQUIRE(a + b == c);
    }

    {
        BivariateRingPolynomial<int> c{3, 2};
        c.coeffs = {12, 6, 27, 2, 2, 2};
        REQUIRE(a - b == c);
    }

    {
        BivariateRingPolynomial<int> c{3, 2};
        c.coeffs = {-12, 5, -176, -129, -231, -7};
        REQUIRE(a * b == c);
    }
}
