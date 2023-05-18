// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include <cassert>
#include <iostream>
#include <vector>

// Elements of the ring Z_q[Z] / (Z^n + 1)
template <typename T>
struct UnivariateRingPolynomial {
    std::size_t n{};
    std::vector<T> coeffs;

    // Make this default constructable
    explicit UnivariateRingPolynomial() = default;

    explicit UnivariateRingPolynomial(std::size_t n_) : n(n_), coeffs(n) {}

    friend std::ostream& operator<<(std::ostream& os, const UnivariateRingPolynomial& obj) {
        bool first = true;
        for (std::size_t i = 0; i < obj.n; ++i) {
            if (!first) {
                os << " + ";
            }
            os << obj.coeffs[i];
            if (i != 0) {
                os << "z";
                if (i != 1) {
                    os << "^" << i;
                }
            }
            first = false;
        }
        return os;
    }

    UnivariateRingPolynomial& operator+=(const UnivariateRingPolynomial& other) {
        assert(n == other.n);
        for (std::size_t i = 0; i < n; ++i) {
            coeffs[i] += other.coeffs[i];
        }
        return *this;
    }

    friend UnivariateRingPolynomial operator+(UnivariateRingPolynomial lhs,
                                              const UnivariateRingPolynomial& rhs) {
        lhs += rhs;
        return lhs;
    }

    UnivariateRingPolynomial& operator+=(const T& other) {
        coeffs[0] += other;
        return *this;
    }

    friend UnivariateRingPolynomial operator+(UnivariateRingPolynomial lhs, const T& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend UnivariateRingPolynomial operator+(const T& lhs, UnivariateRingPolynomial rhs) {
        rhs += lhs;
        return rhs;
    }

    UnivariateRingPolynomial& operator-=(const UnivariateRingPolynomial& other) {
        assert(n == other.n);
        for (std::size_t i = 0; i < n; ++i) {
            coeffs[i] -= other.coeffs[i];
        }
        return *this;
    }

    friend UnivariateRingPolynomial operator-(UnivariateRingPolynomial lhs,
                                              const UnivariateRingPolynomial& rhs) {
        lhs -= rhs;
        return lhs;
    }

    UnivariateRingPolynomial& operator-=(const T& other) {
        coeffs[0] -= other;
        return *this;
    }

    friend UnivariateRingPolynomial operator-(UnivariateRingPolynomial lhs, const T& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend UnivariateRingPolynomial operator-(const T& lhs, UnivariateRingPolynomial rhs) {
        rhs.coeffs[0] = lhs - rhs.coeffs[0];
        for (std::size_t i = 1; i < rhs.coeffs.size(); ++i) {
            rhs.coeffs[i] = -rhs.coeffs[i];
        }
        return rhs;
    }

    UnivariateRingPolynomial& operator*=(const UnivariateRingPolynomial& other) {
        assert(n == other.n);

        std::vector<T> new_coeffs(n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                if (i + j < n) {
                    new_coeffs[i + j] += coeffs[i] * other.coeffs[j];
                } else {
                    new_coeffs[i + j - n] -= coeffs[i] * other.coeffs[j];
                }
            }
        }
        coeffs.swap(new_coeffs);
        return *this;
    }

    friend UnivariateRingPolynomial operator*(UnivariateRingPolynomial lhs,
                                              const UnivariateRingPolynomial& rhs) {
        lhs *= rhs;
        return lhs;
    }

    UnivariateRingPolynomial& operator*=(const T& other) {
        for (std::size_t i = 0; i < n; ++i) {
            coeffs[i] *= other;
        }
        return *this;
    }

    friend UnivariateRingPolynomial operator*(UnivariateRingPolynomial lhs, const T& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend UnivariateRingPolynomial operator*(const T& lhs, UnivariateRingPolynomial rhs) {
        rhs *= lhs;
        return rhs;
    }

    friend bool operator==(const UnivariateRingPolynomial& lhs,
                           const UnivariateRingPolynomial& rhs) {
        assert(lhs.n == rhs.n);
        for (std::size_t i = 0; i < lhs.n; ++i) {
            if (lhs.coeffs[i] != rhs.coeffs[i]) {
                return false;
            }
        }
        return true;
    }

    friend bool operator!=(const UnivariateRingPolynomial& lhs,
                           const UnivariateRingPolynomial& rhs) {
        return !(lhs == rhs);
    }
};

// Elements of the ring Z_q[Z, Y] / (Z^n + 1, Y^D + 1)
template <typename T>
struct BivariateRingPolynomial {
    std::size_t n{}, D{};
    std::vector<T> coeffs;

    // Make this default constructable
    explicit BivariateRingPolynomial() = default;

    explicit BivariateRingPolynomial(std::size_t n_, std::size_t D_)
        : n(n_), D(D_), coeffs(n * D) {}

    T& Coeff(std::size_t i, std::size_t j) {
        return coeffs[i * D + j];
    }

    const T& Coeff(std::size_t i, std::size_t j) const {
        return coeffs[i * D + j];
    }

    friend std::ostream& operator<<(std::ostream& os, const BivariateRingPolynomial& obj) {
        bool first = true;
        for (std::size_t i = 0; i < obj.n; ++i) {
            for (std::size_t j = 0; j < obj.D; ++j) {
                if (!first) {
                    os << " + ";
                }
                os << obj.Coeff(i, j);
                if (i != 0) {
                    os << "z";
                    if (i != 1) {
                        os << "^" << i;
                    }
                }
                if (j != 0) {
                    os << "y";
                    if (j != 1) {
                        os << "^" << j;
                    }
                }
                first = false;
            }
        }
        return os;
    }

    BivariateRingPolynomial& operator+=(const BivariateRingPolynomial& other) {
        assert(n == other.n && D == other.D);
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            coeffs[i] += other.coeffs[i];
        }
        return *this;
    }

    friend BivariateRingPolynomial operator+(BivariateRingPolynomial lhs,
                                             const BivariateRingPolynomial& rhs) {
        lhs += rhs;
        return lhs;
    }

    BivariateRingPolynomial& operator-=(const BivariateRingPolynomial& other) {
        assert(n == other.n && D == other.D);
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            coeffs[i] -= other.coeffs[i];
        }
        return *this;
    }

    friend BivariateRingPolynomial operator-(BivariateRingPolynomial lhs,
                                             const BivariateRingPolynomial& rhs) {
        lhs -= rhs;
        return lhs;
    }

    BivariateRingPolynomial& operator*=(const BivariateRingPolynomial& other) {
        assert(n == other.n && D == other.D);

        std::vector<T> new_coeffs(n * D);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < D; ++j) {
                for (std::size_t k = 0; k < n; ++k) {
                    for (std::size_t l = 0; l < D; ++l) {
                        const bool flip = (i + k < n) != (j + l < D);
                        if (flip) {
                            new_coeffs[((i + k) % n) * D + (j + l) % D] -=
                                Coeff(i, j) * other.Coeff(k, l);
                        } else {
                            new_coeffs[((i + k) % n) * D + (j + l) % D] +=
                                Coeff(i, j) * other.Coeff(k, l);
                        }
                    }
                }
            }
        }
        coeffs.swap(new_coeffs);
        return *this;
    }

    friend BivariateRingPolynomial operator*(BivariateRingPolynomial lhs,
                                             const BivariateRingPolynomial& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend bool operator==(const BivariateRingPolynomial& lhs, const BivariateRingPolynomial& rhs) {
        assert(lhs.n == rhs.n && lhs.D == rhs.D);
        for (std::size_t i = 0; i < lhs.coeffs.size(); ++i) {
            if (lhs.coeffs[i] != rhs.coeffs[i]) {
                return false;
            }
        }
        return true;
    }

    friend bool operator!=(const BivariateRingPolynomial& lhs, const BivariateRingPolynomial& rhs) {
        return !(lhs == rhs);
    }
};
