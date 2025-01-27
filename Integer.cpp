#include <atomic>
#include <mutex>
#include <thread>

#include "Integer.h"

constexpr std::array<unsigned int, 3> const primes{3, 5, 7};

Integer::Integer() : bits_{0}
{
}

Integer::Integer(std::vector<uintmax_t> const& bits, bool isPositive) : isPositive_{isPositive}, bits_{bits}
{
    if (bits_.empty())
        bits_.emplace_back(0);

    adjust();
}

CONSTEXPR std::vector<uintmax_t> const& Integer::bits() const noexcept
{
    return bits_;
}

CONSTEXPR bool Integer::isPositive() const noexcept
{
    return isPositive_;
}

Integer Integer::operator-() const
{
    auto x(*this);

    x.isPositive_ = !x.isPositive_;

    return x;
}

Integer Integer::operator~() const
{
    auto x(*this);

    x.invert();

    return x;
}

Integer::operator bool() const noexcept
{
    return !!*this;
}

Integer& Integer::operator--()
{
    return *this -= 1;
}

Integer Integer::operator--(int)
{
    auto x(*this);

    operator--();

    return x;
}

Integer& Integer::operator++()
{
    return *this += 1;
}

Integer Integer::operator++(int)
{
    auto x(*this);

    operator++();

    return x;
}

CONSTEXPR Integer::operator char() const noexcept
{
    return cast<char>();
}

CONSTEXPR Integer::operator unsigned char() const noexcept
{
    return cast<unsigned char>();
}

CONSTEXPR Integer::operator short() const noexcept
{
    return cast<short>();
}

CONSTEXPR Integer::operator unsigned short() const noexcept
{
    return cast<unsigned short>();
}

CONSTEXPR Integer::operator int() const noexcept
{
    return cast<int>();
}

CONSTEXPR Integer::operator unsigned int() const noexcept
{
    return cast<unsigned int>();
}

CONSTEXPR Integer::operator long() const noexcept
{
    return cast<long>();
}

CONSTEXPR Integer::operator unsigned long() const noexcept
{
    return cast<unsigned long>();
}

CONSTEXPR Integer::operator long long() const noexcept
{
    return cast<long long>();
}

CONSTEXPR Integer::operator unsigned long long() const noexcept
{
    return cast<unsigned long long>();
}

CONSTEXPR bool Integer::isNan() const noexcept
{
    return isNan_;
}

void Integer::setNan() noexcept
{
    isNan_ = true;
    isInfinity_ = false;
    bits_.clear();
}

CONSTEXPR bool Integer::isInfinity() const noexcept
{
    return isInfinity_;
}

void Integer::setInfinity() noexcept
{
    isNan_ = false;
    isInfinity_ = true;
    bits_.clear();
}

Integer Integer::abs() const
{
    if (isNegative())
        return -*this;

    return *this;
}

Integer Integer::number() const noexcept
{
    Integer number(0);

    if (isNan() || isInfinity())
        return number;

    auto it{bits_.begin()};

    while (!*it && it != bits_.end())
        ++it;

    if (it != bits_.end())
    {
        auto b{*it};

        while (b)
        {
            ++number;
            b >>= 1;
        }

        number += (std::distance(it, bits_.end()) - 1) * sizeof(uintmax_t) * 8;
    }

    return number;
}

Integer Integer::sign() const
{
    if (*this < 0)
        return Integer(-1);

    return Integer(1);
}

size_t Integer::size() const noexcept
{
    return bits_.size();
}

void Integer::adjust()
{
    if (bits_.empty())
        return;

    auto it{bits_.begin()};

    while (!*it && it != bits_.end())
        ++it;

    if (it == bits_.end())
        it = bits_.end() - 1;

    if (it != bits_.begin())
        bits_ = std::vector<uintmax_t>{it, bits_.end()};
}

CONSTEXPR bool Integer::autoAdjust() const noexcept
{
    return autoAdjust_;
}

CONSTEXPR void Integer::setAutoAdjust(bool autoAdjust) noexcept
{
    autoAdjust_ = autoAdjust;
}

char const* Integer::data() const noexcept
{
    return reinterpret_cast<char const*>(bits_.data());
}

char* Integer::data() noexcept
{
    return reinterpret_cast<char*>(bits_.data());
}

size_t Integer::dataSize() const noexcept
{
    return sizeof(uintmax_t) / sizeof(char) * bits_.size();
}

void Integer::setData(char const* data, size_t size) noexcept
{
    std::memcpy(bits_.data(), data, std::min(size, dataSize()));
}
