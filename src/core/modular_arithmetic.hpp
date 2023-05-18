// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include <cassert>
#include <iostream>
#include <memory>

// Modular arithmetic. Optional tag to ensure type safety.
template <typename T, typename Tag = void>
struct ModularInt {
    T value{};
    std::shared_ptr<T> modulus{};

    constexpr explicit ModularInt() = default;

    constexpr explicit ModularInt(const T& value_, const std::shared_ptr<T>& modulus_)
        : value(value_ % *modulus_), modulus(modulus_) {}

    constexpr ModularInt(const ModularInt& other) : value(other.value), modulus(other.modulus) {}

    constexpr ModularInt(ModularInt&& other)
        : value(std::move(other.value)), modulus(std::move(other.modulus)) {}

    ModularInt& operator=(const ModularInt& other) {
        value = other.value;
        modulus = other.modulus;
        return *this;
    }

    ModularInt& operator=(ModularInt&& other) {
        value = std::move(other.value);
        modulus = std::move(other.modulus);
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const ModularInt& obj) {
        if (!obj.modulus) {
            os << "0";
            return os;
        }
        os << obj.GetRepresentative();
        return os;
    }

    constexpr ModularInt& operator+=(const ModularInt& other) {
        if (!modulus) {
            modulus = other.modulus;
        }
        if (!other.modulus) {
            return *this;
        }

        assert(*modulus == *other.modulus);
        value += other.value;
        value %= *modulus;
        return *this;
    }

    friend ModularInt operator+(ModularInt lhs, const ModularInt& rhs) {
        lhs += rhs;
        return lhs;
    }

    constexpr ModularInt& operator-=(const ModularInt& other) {
        if (!modulus) {
            modulus = other.modulus;
        }
        if (!other.modulus) {
            return *this;
        }

        assert(*modulus == *other.modulus);
        value -= other.value;
        value %= *modulus;
        return *this;
    }
    friend ModularInt operator-(ModularInt lhs, const ModularInt& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend ModularInt operator-(const ModularInt& val) {
        return ModularInt{-val.value, val.modulus};
    }

    // Can probably be optimized
    constexpr ModularInt& operator*=(const ModularInt& other) {
        if (!modulus) {
            modulus = other.modulus;
        }
        if (!other.modulus) {
            value = 0;
            return *this;
        }

        assert(*modulus == *other.modulus);
        value *= other.value;
        value %= *modulus;
        return *this;
    }

    friend ModularInt operator*(ModularInt lhs, const ModularInt& rhs) {
        lhs *= rhs;
        return lhs;
    }

    friend bool operator==(const ModularInt& lhs, const ModularInt& rhs) {
        if ((!lhs.modulus && rhs.value == 0) || (!rhs.modulus && lhs.value == 0)) {
            return true;
        } else if (!lhs.modulus || !rhs.modulus) {
            return false;
        }

        assert(*lhs.modulus == *rhs.modulus);
        return lhs.value == rhs.value || lhs.value == rhs.value + *lhs.modulus ||
               lhs.value == rhs.value - *lhs.modulus;
    }

    friend bool operator!=(const ModularInt& lhs, const ModularInt& rhs) {
        return !(lhs == rhs);
    }

    T GetRepresentative() const {
        assert(modulus);

        auto rep = value; // (-M, M)
        // To (-M/2, M/2]
        if (rep > *modulus / 2) {
            rep -= *modulus;
        } else if (rep <= -*modulus / 2) {
            rep += *modulus;
        }
        return rep;
    }

    template <typename NewTag = Tag>
    ModularInt<T, NewTag> Reinterpret(const std::shared_ptr<T>& new_modulus) const {
        return ModularInt<T, NewTag>{GetRepresentative() % *new_modulus, new_modulus};
    }
};
