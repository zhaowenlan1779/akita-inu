// Copyright 2023 Pengfei Zhu
// Licensed under MIT.

#pragma once

#include "core/big_int.h"
#include "core/modular_arithmetic.hpp"
#include "core/ring_polynomial.hpp"

// General typedefs
struct qIntTag {};
using qInt = ModularInt<BigInt, qIntTag>;
struct dIntTag {};
using dInt = ModularInt<int64_t, dIntTag, true>;

using Poly = UnivariatePolynomial<dInt>;
using MultiPoly = MultivariatePolynomial<dInt>;
