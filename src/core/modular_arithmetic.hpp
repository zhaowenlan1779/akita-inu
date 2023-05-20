// Copyright 2023 Pengfei Zhu
// Licensed under GPLv2 or any later version.

#pragma once

#include <cassert>
#include <iostream>
#include <memory>
#include <type_traits>
#include <boost/integer/extended_euclidean.hpp>

// Modular arithmetic. Optional tag to ensure type safety.
template <typename T, typename Tag = void, bool EnableDivision = false>
struct ModularInt {
    using ModulusType = std::conditional_t<std::is_scalar_v<T>, T, std::shared_ptr<T>>;
    ModulusType modulus{};

    T value{};

    constexpr bool HasModulus() const {
        return static_cast<bool>(modulus);
    }

    constexpr const T& Modulus() const {
        if constexpr (std::is_scalar_v<T>) {
            return modulus;
        } else {
            return *modulus;
        }
    }

    constexpr T& Modulus() {
        if constexpr (std::is_scalar_v<T>) {
            return modulus;
        } else {
            return *modulus;
        }
    }

    constexpr explicit ModularInt() = default;

    constexpr ModularInt(const T& value_, const ModulusType& modulus_)
        : modulus(modulus_), value(value_ % Modulus()) {}

    constexpr ModularInt(const ModularInt& other) : modulus(other.modulus), value(other.value) {}

    constexpr ModularInt(ModularInt&& other)
        : modulus(std::move(other.modulus)), value(std::move(other.value)) {}

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
        if (!obj.HasModulus()) {
            os << "0";
            return os;
        }
        os << obj.GetRepresentative();
        return os;
    }

    constexpr ModularInt& operator+=(const ModularInt& other) {
        if (!HasModulus()) {
            modulus = other.modulus;
        }
        if (!other.HasModulus()) {
            return *this;
        }

        assert(Modulus() == other.Modulus());
        value += other.value;
        value %= Modulus();
        return *this;
    }

    friend ModularInt operator+(ModularInt lhs, const ModularInt& rhs) {
        lhs += rhs;
        return lhs;
    }

    constexpr ModularInt& operator-=(const ModularInt& other) {
        if (!HasModulus()) {
            modulus = other.modulus;
        }
        if (!other.HasModulus()) {
            return *this;
        }

        assert(Modulus() == other.Modulus());
        value -= other.value;
        value %= Modulus();
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
        if (!HasModulus()) {
            modulus = other.modulus;
        }
        if (!other.HasModulus()) {
            value = 0;
            return *this;
        }

        assert(Modulus() == other.Modulus());
        // TODO: This may get out of range with scalars
        value *= other.value;
        value %= Modulus();
        return *this;
    }

    friend ModularInt operator*(ModularInt lhs, const ModularInt& rhs) {
        lhs *= rhs;
        return lhs;
    }

    constexpr ModularInt& operator/=(const ModularInt& other) {
        if constexpr (!EnableDivision) {
            char x[0] = "Division was not enabled";
        }
        if (!HasModulus()) {
            modulus = other.modulus;
        }
        if (!other.HasModulus()) {
            assert(false);
            return *this;
        }

        assert(Modulus() == other.Modulus());

        // Find inverse of other.value
        value *= boost::integer::extended_euclidean(other.value + Modulus(), Modulus()).x;
        value %= Modulus();
        return *this;
    }

    friend ModularInt operator/(ModularInt lhs, const ModularInt& rhs) {
        lhs /= rhs;
        return lhs;
    }

    friend bool operator==(const ModularInt& lhs, const ModularInt& rhs) {
        if ((!lhs.HasModulus() && rhs.value == 0) || (!rhs.HasModulus() && lhs.value == 0)) {
            return true;
        } else if (!lhs.HasModulus() || !rhs.HasModulus()) {
            return false;
        }

        assert(lhs.Modulus() == rhs.Modulus());
        return lhs.value == rhs.value || lhs.value == rhs.value + lhs.Modulus() ||
               lhs.value == rhs.value - lhs.Modulus();
    }

    friend bool operator!=(const ModularInt& lhs, const ModularInt& rhs) {
        return !(lhs == rhs);
    }

    T GetRepresentative() const {
        assert(HasModulus());

        auto rep = value; // (-M, M)
        // To (-M/2, M/2]
        if (rep > Modulus() / 2) {
            rep -= Modulus();
        } else if (rep <= -Modulus() / 2) {
            rep += Modulus();
        }
        return rep;
    }

    template <typename NewTag = Tag>
    ModularInt<T, NewTag> Reinterpret(const ModulusType& new_modulus) const {
        if constexpr (std::is_scalar_v<T>) {
            return ModularInt<T, NewTag>{GetRepresentative() % new_modulus, new_modulus};
        } else {
            return ModularInt<T, NewTag>{GetRepresentative() % *new_modulus, new_modulus};
        }
    }

    explicit constexpr operator bool() const {
        return HasModulus() && value != 0;
    }
};
