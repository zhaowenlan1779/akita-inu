// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include <cassert>
#include <iostream>

// Modular arithmetic.
template <typename T>
struct ModularInt {
    T value;
    const T& modulus;

    constexpr explicit ModularInt(const T& value_, const T& modulus_)
        : value(value_ % modulus_), modulus(modulus_) {}

    friend std::ostream& operator<<(std::ostream& os, const ModularInt& obj) {
        os << obj.GetRepresentative();
        return os;
    }

    constexpr ModularInt& operator+=(const ModularInt& other) {
        assert(modulus == other.modulus);
        value += other.value;
        value %= modulus;
        return *this;
    }

    constexpr ModularInt& operator+=(const T& other) {
        value += other;
        value %= modulus;
        return *this;
    }

    friend ModularInt operator+(ModularInt lhs, const ModularInt& rhs) {
        lhs += rhs;
        return lhs;
    }

    friend ModularInt operator+(const T& lhs, ModularInt rhs) {
        rhs += lhs;
        return rhs;
    }

    friend ModularInt operator+(ModularInt lhs, const T& rhs) {
        lhs += rhs;
        return lhs;
    }

    constexpr ModularInt& operator-=(const ModularInt& other) {
        assert(modulus == other.modulus);
        value -= other.value;
        value %= modulus;
        return *this;
    }

    constexpr ModularInt& operator-=(const T& other) {
        value -= other;
        value %= modulus;
        return *this;
    }

    friend ModularInt operator-(ModularInt lhs, const ModularInt& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend ModularInt operator-(const T& lhs, ModularInt rhs) {
        rhs.value = lhs - rhs.value;
        rhs.value %= rhs.modulus;
        return rhs;
    }

    friend ModularInt operator-(ModularInt lhs, const T& rhs) {
        lhs -= rhs;
        return lhs;
    }

    // Can probably be optimized
    constexpr ModularInt& operator*=(const ModularInt& other) {
        assert(modulus == other.modulus);
        value *= other.value;
        value %= modulus;
        return *this;
    }

    constexpr ModularInt& operator*=(const T& other) {
        value *= other;
        value %= modulus;
        return *this;
    }

    friend ModularInt operator*(ModularInt lhs, const ModularInt& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend ModularInt operator*(const T& lhs, ModularInt rhs) {
        rhs *= lhs;
        return rhs;
    }

    friend ModularInt operator*(ModularInt lhs, const T& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend bool operator==(const ModularInt& lhs, const ModularInt& rhs) {
        assert(lhs.modulus == rhs.modulus);
        return lhs.value == rhs.value || lhs.value == rhs.value + lhs.modulus ||
               lhs.value == rhs.value - lhs.modulus;
    }

    friend bool operator!=(const ModularInt& lhs, const ModularInt& rhs) {
        return !(lhs == rhs);
    }

    T GetRepresentative() const {
        // As described in the paper, the representative is in (-M/2, M/2].
        // First convert to (-M, 0], then to (-M/2, M/2] by +M/2.
        // Note that the original value may be in (-M, M)
        return (value - modulus) % modulus + modulus / 2;
    }

    ModularInt Reinterpret(const T& new_modulus) const {
        return ModularInt{GetRepresentative() % new_modulus, new_modulus};
    }
};
