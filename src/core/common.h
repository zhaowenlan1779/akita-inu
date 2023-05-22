// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include "core/modular_arithmetic.hpp"
#include "core/ring_polynomial.hpp"

// General typedefs
struct dIntTag {};
using dInt = ModularInt<int64_t, dIntTag, true>;

using Poly = UnivariatePolynomial<dInt>;
using MultiPoly = MultivariatePolynomial<dInt>;
