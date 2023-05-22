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
    using ValueType = T;
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

    constexpr ModularInt(const T& value_, const ModulusType& modulus_) noexcept
        : modulus(modulus_), value(value_ % Modulus()) {}

    constexpr ModularInt(const T& value_) noexcept : value(value_) {}

    constexpr ModularInt(T&& value_) noexcept : value(value_) {}

    constexpr ModularInt(const ModularInt& other) noexcept
        : modulus(other.modulus), value(other.value) {}

    constexpr ModularInt(ModularInt&& other) noexcept
        : modulus(std::move(other.modulus)), value(std::move(other.value)) {}

    ModularInt& operator=(const ModularInt& other) noexcept {
        value = other.value;
        modulus = other.modulus;
        return *this;
    }

    ModularInt& operator=(ModularInt&& other) noexcept {
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
        if (HasModulus() && other.HasModulus()) {
            assert(Modulus() == other.Modulus());
        } else if (!HasModulus()) {
            modulus = other.modulus;
        }
        value += other.value;
        if (HasModulus()) {
            value %= Modulus();
        }
        return *this;
    }

    friend ModularInt operator+(ModularInt lhs, const ModularInt& rhs) {
        lhs += rhs;
        return lhs;
    }

    constexpr ModularInt& operator-=(const ModularInt& other) {
        if (HasModulus() && other.HasModulus()) {
            assert(Modulus() == other.Modulus());
        } else if (!HasModulus()) {
            modulus = other.modulus;
        }
        value -= other.value;
        if (HasModulus()) {
            value %= Modulus();
        }
        return *this;
    }

    friend ModularInt operator-(ModularInt lhs, const ModularInt& rhs) {
        lhs -= rhs;
        return lhs;
    }

    friend ModularInt operator-(const ModularInt& val) {
        ModularInt new_val{val};
        new_val.value = -new_val.value;
        return new_val;
    }

    // Can probably be optimized
    constexpr ModularInt& operator*=(const ModularInt& other) {
        if (HasModulus() && other.HasModulus()) {
            assert(Modulus() == other.Modulus());
        } else if (!HasModulus()) {
            modulus = other.modulus;
        }
        // TODO: This may get out of range with scalars
        value *= other.value;
        if (HasModulus()) {
            value %= Modulus();
        }
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

        if (HasModulus() && other.HasModulus()) {
            assert(Modulus() == other.Modulus());
        } else if (!HasModulus()) {
            modulus = other.modulus;
        }
        assert(HasModulus());

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
        if (!lhs.HasModulus() && !rhs.HasModulus()) {
            return lhs.value == rhs.value;
        }

        assert(lhs.HasModulus() && rhs.HasModulus());
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
        } else if (rep <= -(Modulus() + 1) / 2) {
            rep += Modulus();
        }
        return rep;
    }

    template <typename U = ModularInt>
    U Reinterpret(const typename U::ModulusType& new_modulus) const {
        using NewModulusType = U::ModulusType;
        using NewValueType = U::ValueType;
        using IntermediateType =
            std::conditional_t<std::conjunction_v<std::is_scalar<T>, std::is_scalar<NewValueType>>,
                               std::common_type_t<T, NewValueType>,
                               std::conditional_t<std::is_scalar_v<T>, NewValueType, T>>;

        if constexpr (std::is_scalar_v<NewModulusType>) {
            return U{
                NewValueType{IntermediateType{GetRepresentative()} % IntermediateType{new_modulus}},
                new_modulus};
        } else {
            return U{NewValueType{IntermediateType{GetRepresentative()} %
                                  IntermediateType{*new_modulus}},
                     new_modulus};
        }
    }

    explicit constexpr operator bool() const {
        return HasModulus() && value != 0;
    }
};
