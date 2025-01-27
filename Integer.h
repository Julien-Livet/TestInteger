#ifndef INTEGER_H
#define INTEGER_H

/*!
 *  \file Integer.h
 *  \brief Provide a class to manage large integer numbers
 *  \author Julien LIVET
 *  \version 1.0
 *  \date 28/12/2024
 */

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstring>
#include <functional>
#include <iomanip>
#include <random>
#include <sstream>
#include <vector>

#if __cplusplus >= 201703L
#define CONSTEXPR constexpr
#else
#define CONSTEXPR
#endif

#if __cplusplus >= 202106L
#include <format>
#endif

#ifdef WITH_GMP
#include <gmpxx.h>
#endif

template <typename E>
class IntegerExpression
{
    public:
        static constexpr bool is_leaf = false;

        CONSTEXPR bool isPositive() const
        {
            return static_cast<E const&>(*this).isPositive();
        }

        CONSTEXPR bool isNegative() const
        {
            return !isPositive();
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            return static_cast<E const&>(*this).bits();
        }

        CONSTEXPR bool isNan() const
        {
            return static_cast<E const&>(*this).isNan();
        }

        CONSTEXPR bool isInfinity() const
        {
            return static_cast<E const&>(*this).isInfinity();
        }

        CONSTEXPR operator bool() const noexcept
        {
            return !!*this;
        }

        CONSTEXPR bool operator!() const noexcept
        {
            for (auto const& b : bits())
            {
                if (b)
                    return false;
            }

            return true;
        }

        template <typename T>
        CONSTEXPR bool operator==(IntegerExpression<T> const& other) const
        {
            if (isNan() && other.isNan())
                return true;
            else if (isNan() || other.isNan())
                return false;
            else if (isInfinity() && other.isInfinity())
                return (isPositive() == other.isPositive());
            else if (isInfinity() || other.isInfinity())
                return false;

            if (bits().size() != other.bits().size())
            {
                if (bits().size() > other.bits().size())
                {
                    for (size_t i{0}; i < bits().size() - other.bits().size(); ++i)
                    {
                        if (bits()[i])
                            return false;
                    }
                }
                else
                {
                    for (size_t i{0}; i < other.bits().size() - bits().size(); ++i)
                    {
                        if (other.bits()[i])
                            return false;
                    }
                }
            }

            bool zero{true};

            auto it1{bits().rbegin()};
            auto it2{other.bits().rbegin()};

            for (size_t i{0}; i < std::min(bits().size(), other.bits().size()); ++i)
            {
                if (*it1 != *it2)
                    return false;

                if (*it1)
                    zero = false;

                ++it1;
                ++it2;
            }

            if (isPositive() != other.isPositive() && !zero)
                return false;

            return true;
        }

        template <typename T>
        CONSTEXPR bool operator!=(IntegerExpression<T> const& other) const
        {
            return !operator==(other);
        }

        template <typename T>
        CONSTEXPR bool operator<=(IntegerExpression<T> const& other) const
        {
            return operator<(other) || operator==(other);
        }

        template <typename T>
        CONSTEXPR bool operator>=(IntegerExpression<T> const& other) const
        {
            return operator>(other) || operator==(other);
        }

        template <typename T>
        CONSTEXPR bool operator<(IntegerExpression<T> const& other) const
        {
            if (isPositive() && !other.isPositive())
                return false;
            else if (isNan() || other.isNan())
                return false;
            else if (other.isInfinity())
            {
                if (other.isPositive())
                {
                    if (isInfinity() && isPositive())
                        return false;
                    else
                        return true;
                }
                else
                    return false;
            }

            std::vector<uintmax_t> a(std::max(bits().size(), other.bits().size()), 0);
            std::vector<uintmax_t> b{a};

            std::copy(bits().rbegin(), bits().rend(), a.rbegin());
            std::copy(other.bits().rbegin(), other.bits().rend(), b.rbegin());

            auto const less{a < b};

            return isPositive() ? less : !less;
        }

        template <typename T>
        CONSTEXPR bool operator>(IntegerExpression<T> const& other) const
        {
            if (!isPositive() && other.isPositive())
                return false;
            else if (isNan() || other.isNan())
                return false;
            else if (other.isInfinity())
            {
                if (other.isNegative())
                {
                    if (isInfinity() && isNegative())
                        return false;
                    else
                        return true;
                }
                else
                    return false;
            }

            std::vector<uintmax_t> a(std::max(bits().size(), other.bits().size()), 0);
            std::vector<uintmax_t> b{a};

            std::copy(bits().rbegin(), bits().rend(), a.rbegin());
            std::copy(other.bits().rbegin(), other.bits().rend(), b.rbegin());

            auto const great{a > b};

            return isPositive() ? great : !great;
        }

        template <typename S>
        CONSTEXPR S cast() const
        {
            return static_cast<E const&>(*this).template cast<S>();
        }

        template <typename S>
        CONSTEXPR bool fits() const
        {
            return static_cast<E const&>(*this).template fits<S>();
        }
};

class Integer : public IntegerExpression<Integer>
{
    public:
        static constexpr bool is_leaf = true;

        Integer();

        template <typename E>
        Integer(IntegerExpression<E> const& expr)
        {
            isPositive_ = expr.isPositive();
            bits_ = expr.bits();
            isNan_ = expr.isNan();
            isInfinity_ = expr.isInfinity();
        }

        template <typename S, std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
        CONSTEXPR Integer(S n) : isPositive_{n >= 0}
        {
            bits_.reserve(1);

            if (n < 0)
                n = -n;

            bits_.emplace_back(n);

            adjust();
        }

        Integer(std::vector<uintmax_t> const& bits, bool isPositive = true);
        CONSTEXPR std::vector<uintmax_t> const& bits() const noexcept;
        void invert() noexcept;
        CONSTEXPR bool isPositive() const noexcept;

        template <typename S>
        CONSTEXPR Integer& operator*=(S const& other)
        {
            *this = *this * other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator+=(S const& other)
        {
            *this = *this + other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator-=(S const& other)
        {
            *this = *this - other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator/=(S const& other)
        {
            *this = *this / other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator%=(S const& other)
        {
            *this = *this % other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator<<=(S const& other)
        {
            *this = *this << other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator>>=(S const& other)
        {
            *this = *this >> other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator&=(S const& other)
        {
            *this = *this & other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator|=(S const& other)
        {
            *this = *this | other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator^=(S const& other)
        {
            *this = *this ^ other;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR bool operator>=(S const& other) const
        {
            return operator>(Integer(other)) || operator==(Integer(other));
        }

        template <typename S>
        CONSTEXPR bool operator>(S const& other) const
        {
            return operator>(Integer(other));
        }

        template <typename S>
        CONSTEXPR bool operator<=(S const& other) const
        {
            return operator<(Integer(other)) || operator==(Integer(other));
        }

        template <typename S>
        CONSTEXPR bool operator<(S const& other) const
        {
            return operator<(Integer(other));
        }

        template <typename S>
        CONSTEXPR bool operator==(S const& other) const
        {
            return *this == Integer(other);
        }

        template <typename S>
        CONSTEXPR bool operator!=(S const& other) const
        {
            return *this != Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator=(S const& other)
        {
            return *this = Integer(other);
        }

        template <typename S>
        CONSTEXPR S cast() const
        {
            if (bits_.empty())
                return S{0};

            auto n{static_cast<S>(bits_.back())};

            if (!isPositive_)
                n = -n;

            return n;
        }

        Integer operator-() const;
        Integer operator~() const;
        explicit operator bool() const noexcept;
        Integer& operator--();
        Integer operator--(int);
        Integer& operator++();
        Integer operator++(int);
        CONSTEXPR operator char() const noexcept;
        CONSTEXPR operator unsigned char() const noexcept;
        CONSTEXPR operator short() const noexcept;
        CONSTEXPR operator unsigned short() const noexcept;
        CONSTEXPR operator int() const noexcept;
        CONSTEXPR operator unsigned int() const noexcept;
        CONSTEXPR operator long() const noexcept;
        CONSTEXPR operator unsigned long() const noexcept;
        CONSTEXPR operator long long() const noexcept;
        CONSTEXPR operator unsigned long long() const noexcept;
        CONSTEXPR bool isNan() const noexcept;
        void setNan() noexcept;
        CONSTEXPR bool isInfinity() const noexcept;
        void setInfinity() noexcept;
        Integer abs() const;

        static Integer nan()
        {
            static Integer n;
            n.setNan();

            return n;
        }

        static Integer infinity()
        {
            static Integer n;
            n.setInfinity();

            return n;
        }

        Integer number() const noexcept;

        template <typename S>
        CONSTEXPR bool fits() const
        {
            return (*this == this->cast<S>());
        }

        Integer sign() const;
        size_t size() const noexcept;
        void adjust();
        CONSTEXPR bool autoAdjust() const noexcept;
        CONSTEXPR void setAutoAdjust(bool autoAdjust) noexcept;
        char const* data() const noexcept;
        char* data() noexcept;
        size_t dataSize() const noexcept;
        void setData(char const* data, size_t size) noexcept;

    private:
        bool isPositive_{true};
        std::vector<uintmax_t> bits_;
        bool isNan_{false};
        bool isInfinity_{false};
        bool autoAdjust_{true};
};

inline Integer number(Integer const& n) noexcept
{
    return n.number();
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator>(S const& lhs, IntegerExpression<E> const& rhs) noexcept
{
    return rhs.operator<(Integer(lhs));
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator>=(S const& lhs, IntegerExpression<E> const& rhs) noexcept
{
    return rhs.operator<=(Integer(lhs));
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator<(S const& lhs, IntegerExpression<E> const& rhs) noexcept
{
    return rhs.operator>(Integer(lhs));
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator<=(S const& lhs, IntegerExpression<E> const& rhs) noexcept
{
    return rhs.operator>=(Integer(lhs));
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator==(S const& lhs, IntegerExpression<E> const& rhs) noexcept
{
    return rhs.operator==(Integer(lhs));
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator!=(S const& lhs, IntegerExpression<E> const& rhs) noexcept
{
    return rhs.operator!=(Integer(lhs));
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator==(IntegerExpression<E> const& lhs, S const& rhs) noexcept
{
    return lhs.operator==(Integer(rhs));
}

template <typename E, typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline bool operator!=(IntegerExpression<E> const& lhs, S const& rhs) noexcept
{
    return lhs.operator!=(Integer(rhs));
}

inline Integer abs(Integer const& n)
{
    return n.abs();
}

template <typename E1, typename E2>
class IntegerAnd : public IntegerExpression<IntegerAnd<E1, E2> >
{
    public:
        CONSTEXPR IntegerAnd(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            return u_.isPositive() & v_.isPositive();
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            std::vector<uintmax_t> bits(std::max(u_.bits().size(), v_.bits().size())
                                        - std::min(u_.bits().size(), v_.bits().size()), 0);

            if (u_.bits().size() > v_.bits().size())
                bits.insert(bits.end(), v_.bits().begin(), v_.bits().end());
            else
                bits.insert(bits.end(), u_.bits().begin(), u_.bits().end());

            std::vector<uintmax_t> const other(u_.bits().size() > v_.bits().size() ? u_.bits() : v_.bits());

            auto threadFunc
                {
                    [&bits, other] (size_t start, size_t end) -> void
                    {
                        for (size_t i{start}; i < end; ++i)
                            *(bits.rbegin() + i) &= *(other.rbegin() + i);
                    }
                };

            size_t const numThreads{std::thread::hardware_concurrency()};
            size_t const chunkSize{std::min(u_.bits().size(), v_.bits().size()) / numThreads};
            std::vector<std::thread> threads;

            for (size_t i{0}; i < numThreads; ++i)
            {
                size_t const start{i * chunkSize};
                size_t const end{(i == numThreads - 1) ? std::min(u_.bits().size(), v_.bits().size()) : (i + 1) * chunkSize};
                threads.emplace_back(threadFunc, start, end);
            }

            for (auto& t : threads)
                t.join();

            return bits;
        }

        CONSTEXPR bool isNan() const
        {
            return (u_.isNan() || v_.isNan());
        }

        CONSTEXPR bool isInfinity() const
        {
            return (u_.isInfinity() || v_.isInfinity());
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerAnd<E1, E2> operator&(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerAnd<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
class IntegerOr : public IntegerExpression<IntegerOr<E1, E2> >
{
    public:
        CONSTEXPR IntegerOr(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            return u_.isPositive() | v_.isPositive();
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            std::vector<uintmax_t> bits(u_.bits().size() > v_.bits().size() ? u_.bits() : v_.bits());
            std::vector<uintmax_t> const other(u_.bits().size() > v_.bits().size() ? v_.bits() : u_.bits());

            auto threadFunc
                {
                    [&bits, other] (size_t start, size_t end) -> void
                    {
                        for (size_t i{start}; i < end; ++i)
                            *(bits.rbegin() + i) |= *(other.rbegin() + i);
                    }
                };

            size_t const numThreads{std::thread::hardware_concurrency()};
            size_t const chunkSize{other.size() / numThreads};
            std::vector<std::thread> threads;

            for (size_t i{0}; i < numThreads; ++i)
            {
                size_t const start{i * chunkSize};
                size_t const end{(i == numThreads - 1) ? other.size() : (i + 1) * chunkSize};
                threads.emplace_back(threadFunc, start, end);
            }

            for (auto& t : threads)
                t.join();

            return bits;
        }

        CONSTEXPR bool isNan() const
        {
            return (u_.isNan() || v_.isNan());
        }

        CONSTEXPR bool isInfinity() const
        {
            return (u_.isInfinity() || v_.isInfinity());
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerOr<E1, E2> operator|(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerOr<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
class IntegerXor : public IntegerExpression<IntegerXor<E1, E2> >
{
    public:
        CONSTEXPR IntegerXor(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            return u_.isPositive() | v_.isPositive();
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            std::vector<uintmax_t> bits(u_.bits().size() > v_.bits().size() ? u_.bits() : v_.bits());
            std::vector<uintmax_t> const other(u_.bits().size() > v_.bits().size() ? v_.bits() : u_.bits());

            auto threadFunc1
                {
                    [&bits] (size_t start, size_t end) -> void
                    {
                        for (size_t i{start}; i < end; ++i)
                            *(bits.rbegin() + i) ^= 0;
                    }
                };

            size_t const numThreads{std::thread::hardware_concurrency()};
            size_t const chunkSize{other.size() / numThreads};
            std::vector<std::thread> threads;

            for (size_t i{0}; i < numThreads; ++i)
            {
                size_t const start{i * chunkSize};
                size_t const end{(i == numThreads - 1) ? std::max(bits.size(), other.size())
                                                             - std::min(bits.size(), other.size()) : (i + 1) * chunkSize};
                threads.emplace_back(threadFunc1, start, end);
            }

            for (auto& t : threads)
                t.join();

            auto threadFunc2
                {
                    [&bits, other] (size_t start, size_t end) -> void
                    {
                        for (size_t i{start}; i < end; ++i)
                            *(bits.rbegin() + i) ^= *(other.rbegin() + i);
                    }
                };

            threads.clear();

            for (size_t i{0}; i < numThreads; ++i)
            {
                size_t const start{i * chunkSize};
                size_t const end{(i == numThreads - 1) ? other.size() : (i + 1) * chunkSize};
                threads.emplace_back(threadFunc2, start, end);
            }

            for (auto& t : threads)
                t.join();

            return bits;
        }

        CONSTEXPR bool isNan() const
        {
            return (u_.isNan() || v_.isNan());
        }

        CONSTEXPR bool isInfinity() const
        {
            return (u_.isInfinity() || v_.isInfinity());
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerXor<E1, E2> operator^(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerXor<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator&(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerAnd<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator&(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerAnd<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator|(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerOr<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator|(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerOr<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator^(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerXor<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator^(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerXor<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E1, typename E2>
class IntegerAdd : public IntegerExpression<IntegerAdd<E1, E2> >
{
    public:
        CONSTEXPR IntegerAdd(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            if (u_.isPositive() && v_.isPositive())
                return true;
            else if (!u_.isPositive() && !v_.isPositive())
                return false;
            else if (u_.isPositive() && !v_.isPositive())
                return (u_ > v_);
            else //if (!u_.isPositive() && v_.isPositive())
                return (u_ < v_);
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            if (u_.isNan() || u_.isInfinity() || v_.isNan() || v_.isInfinity())
                return std::vector<uintmax_t>{};
            else
            {
                if ((u_.isPositive() && v_.isPositive())
                    || (u_.isNegative() && v_.isNegative()))
                {
                    uintmax_t carry{0};
                    auto const a{u_.bits()};
                    auto const b{v_.bits()};
                    size_t const n{std::max(a.size(), b.size())};
                    std::vector<uintmax_t> result;
                    result.reserve(n);

                    for (size_t i{0}; i < n; ++i)
                    {
                        auto const bit_a{(i < a.size()) ? a[a.size() - 1 - i] : 0};
                        auto const bit_b{(i < b.size()) ? b[b.size() - 1 - i] : 0};
                        auto const sum{static_cast<uintmax_t>(bit_a + bit_b + carry)};

                        carry = (sum < bit_a || sum < bit_b);

                        result.emplace_back(sum);
                    }

                    if (carry)
                        result.emplace_back(1);

                    std::reverse(result.begin(), result.end());

                    return result;
                }
                else
                {
                    auto a{u_.bits()};
                    auto b{v_.bits()};

                    if (u_.isPositive())
                    {
                        if (u_ < -v_)
                            std::swap(a, b);
                    }
                    else
                    {
                        if (u_ > -v_)
                            std::swap(a, b);
                    }

                    size_t const n{std::max(a.size(), b.size())};
                    std::vector<uintmax_t> result;
                    result.reserve(n);

                    for (size_t i{n - 1}; i <= n - 1; --i)
                    {
                        auto const ia{a.size() - 1 - i};
                        auto const ib{b.size() - 1 - i};

                        auto const bit_a{ia < a.size() ? a[ia] : 0};
                        auto const bit_b{ib < b.size() ? b[ib] : 0};

                        auto bit_result{static_cast<uintmax_t>(bit_a - bit_b)};

                        if (bit_a < bit_b)
                        {
                            for (auto it{result.rbegin()}; it != result.rend(); ++it)
                            {
                                bool const stop{*it > 0};

                                *it -= 1;

                                if (stop)
                                    break;
                            }

                            bit_result = static_cast<uintmax_t>(-1) - bit_b;
                            bit_result += bit_a + 1;
                        }

                        result.emplace_back(bit_result);
                    }

                    return result;
                }
            }
        }

        CONSTEXPR bool isNan() const
        {
            if (u_.isNan() || v_.isNan())
                return true;
            else if ((u_.isInfinity() && u_.isPositive() && v_.isInfinity() && v_.isNegative())
                     || (u_.isInfinity() && u_.isNegative() && v_.isInfinity() && v_.isPositive()))
                return true;
            else
                return false;
        }

        CONSTEXPR bool isInfinity() const
        {
            return (u_.isInfinity() || v_.isInfinity());
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerAdd<E1, E2> operator+(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerAdd<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
class IntegerSub : public IntegerExpression<IntegerSub<E1, E2> >
{
    public:
        CONSTEXPR IntegerSub(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            if (u_.isPositive() && v_.isPositive())
                return (u_ > v_);
            else if (!u_.isPositive() && !v_.isPositive())
                return (u_ < v_);
            else if (u_.isPositive() && !v_.isPositive())
                return true;
            else //if (!u_.isPositive() && v_.isPositive())
                return false;
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            return (u_ + (-v_)).bits();
        }

        CONSTEXPR bool isNan() const
        {
            if (u_.isNan() || v_.isNan())
                return true;
            else if ((u_.isInfinity() && u_.isNegative() && v_.isInfinity() && v_.isNegative())
                     || (u_.isInfinity() && u_.isPositive() && v_.isInfinity() && v_.isNegative()))
                return true;
            else
                return false;
        }

        CONSTEXPR bool isInfinity() const
        {
            return (u_.isInfinity() || v_.isInfinity());
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerSub<E1, E2> operator-(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerSub<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
class IntegerRShift : public IntegerExpression<IntegerRShift<E1, E2> >
{
public:
    CONSTEXPR IntegerRShift(E1 const& u, E2 const& v) : u_(u), v_(v)
    {
    }

    CONSTEXPR bool isPositive() const
    {
        return u_.isPositive();
    }

    CONSTEXPR std::vector<uintmax_t> bits() const
    {
        assert(v_ >= 0);

        if (!u_ || !v_)
            return std::vector<uintmax_t>{0};
        else if (u_.isNan() || v_.isNan() || u_.isInfinity())
            return std::vector<uintmax_t>{};
        else if (v_.isInfinity())
        {
            if (v_ < 0)
                return std::vector<uintmax_t>{};
            else
                return std::vector<uintmax_t>{0};
        }

        auto const s{static_cast<unsigned short>(sizeof(uintmax_t) * 8)};
        auto const n(v_ / s);

        if (u_.bits().size() < n.template cast<uintmax_t>())
            return std::vector<uintmax_t>{0};

        auto bits{u_.bits()};
        bits.resize(bits.size() - n.template cast<uintmax_t>());

        auto const other(v_ - n * s);

        auto const shift{other.template cast<uintmax_t>()};

        if (shift)
        {
            for (auto it{bits.rbegin()}; it != bits.rend(); ++it)
            {
                *it >>= shift;

                if (it != bits.rend() - 1 && (*(it + 1) & ((uintmax_t{1} << shift) - 1)))
                    *it |= (*(it +  1) & ((uintmax_t{1} << shift) - 1)) << (sizeof(uintmax_t) * 8 - shift);
            }
        }

        if (bits.empty())
            bits.emplace_back(0);

        return bits;
    }

    CONSTEXPR bool isNan() const
    {
        if (u_.isNan() || v_.isNan() || u_.isInfinity() || v_.isInfinity())
            return true;
        else
            return false;
    }

    CONSTEXPR bool isInfinity() const
    {
        return u_.isInfinity();
    }

private:
    E1 const u_;
    E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerRShift<E1, E2> operator<<(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerRShift<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
class IntegerLShift : public IntegerExpression<IntegerLShift<E1, E2> >
{
public:
    CONSTEXPR IntegerLShift(E1 const& u, E2 const& v) : u_(u), v_(v)
    {
    }

    CONSTEXPR bool isPositive() const
    {
        return u_.isPositive();
    }

    CONSTEXPR std::vector<uintmax_t> bits() const
    {
        assert(v_ >= 0);

        if (!u_ || !v_)
            return std::vector<uintmax_t>{0};
        else if (u_.isNan() || v_.isNan() || u_.isInfinity() || v_.isInfinity())
            return std::vector<uintmax_t>{};

        auto const s{static_cast<unsigned short>(sizeof(uintmax_t) * 8)};
        auto const n(v_ / s);

        std::vector<uintmax_t> const v(n.template cast<uintmax_t>(), 0);

        auto bits{u_.bits()};
        bits.insert(bits.end(), v.begin(), v.end());

        auto const other(v_ - n * s);

        std::vector<uintmax_t> b(bits.size() + 1, 0);

        std::copy(bits.rbegin(), bits.rend(), b.rbegin());

        bits = b;

        auto const shift{other.template cast<uintmax_t>()};

        if (shift)
        {
            for (auto it{bits.begin() + 1}; it != bits.end(); ++it)
            {
                uintmax_t const s{sizeof(uintmax_t) * 8};

                if ((*it >> (s - shift)))
                    *(it - 1) |= (*it >> (s - shift));

                *it <<= shift;
            }
        }

        return bits;
    }

    CONSTEXPR bool isNan() const
    {
        if (u_.isNan() || v_.isNan() || u_.isInfinity())
            return true;
        else
            return false;
    }

    CONSTEXPR bool isInfinity() const
    {
        return u_.isInfinity();
    }

private:
    E1 const u_;
    E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerLShift<E1, E2> operator>>(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerLShift<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator>>(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerRShift<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator>>(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerRShift<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator<<(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerLShift<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator<<(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerLShift<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E1, typename E2>
class IntegerMul : public IntegerExpression<IntegerMul<E1, E2> >
{
    public:
        CONSTEXPR IntegerMul(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            if (u_.isPositive() && v_.isPositive())
                return true;
            else if (!u_.isPositive() && !v_.isPositive())
                return true;
            else
                return false;
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            if (u_.isNan() || v_.isNan())
                return std::vector<uintmax_t>{};
            else if (!u_ || !v_)
                return std::vector<uintmax_t>{0};
            else if (u_.isInfinity() || v_.isInfinity())
                return std::vector<uintmax_t>{};
            else
            {
                if (u_.isPositive() && v_.isPositive())
                {
                    if (u_.template fits<uintmax_t>() && v_.template fits<uintmax_t>())
                    {
                        auto const a{u_.template cast<uintmax_t>()};
                        auto const b{v_.template cast<uintmax_t>()};
                        auto const ab{a * b};

                        if (ab / b == a)
                            return std::vector<uintmax_t>{ab};
                        else
                        {
                            auto number{[] (uintmax_t n) -> size_t
                                {
                                    size_t number{0};

                                    while (n)
                                    {
                                        ++number;
                                        n >>= 1;
                                    }

                                    return number;
                                }
                            };

                            //Karatsuba algorithm
                            size_t n{std::max(number(a), number(b))};
                            if (n % 2)
                                ++n;
                            size_t const m{n / 2};
                            Integer const x0((static_cast<uintmax_t>(~uintmax_t{0}) >> (sizeof(uintmax_t) * 8 - m)) & a);
                            Integer const x1((((static_cast<uintmax_t>(~uintmax_t{0}) >> (sizeof(uintmax_t) * 8 - m)) << m) & a) >> m);
                            Integer const y0((static_cast<uintmax_t>(~uintmax_t{0}) >> (sizeof(uintmax_t) * 8 - m)) & b);
                            Integer const y1((((static_cast<uintmax_t>(~uintmax_t{0}) >> (sizeof(uintmax_t) * 8 - m)) << m) & b) >> m);

                            auto const z0(IntegerMul<Integer, Integer>(x0, y0));
                            auto const z1(IntegerMul<Integer, Integer>(x1, y0)
                                          + IntegerMul<Integer, Integer>(x0, y1));
                            auto const z2(IntegerMul<Integer, Integer>(x1, y1));

                            //xy = z2 * 2^(2 * m) + z1 * 2^m + z0
                            return (z0 + (z1 << m) + (z2 << 2 * m)).bits();
                        }
                    }
                    else if (!(v_ & 1))
                    {
                        auto r(v_);
                        Integer shift(0);

                        while (!(r & 1))
                        {
                            r >>= 1;
                            ++shift;
                        }

                        return ((u_ << shift) * r).bits();
                    }
                    else
                    {
                        //Karatsuba algorithm
                        //x = x1 * 2^m + x0
                        //y = y1 * 2^m + y0
                        size_t n1{number(u_) / (sizeof(uintmax_t) * 8)};
                        if (number(u_) % (sizeof(uintmax_t) * 8))
                            ++n1;
                        if (n1 % 2)
                            ++n1;
                        size_t n2{number(v_) / (sizeof(uintmax_t) * 8)};
                        if (number(v_) % (sizeof(uintmax_t) * 8))
                            ++n2;
                        if (n2 % 2)
                            ++n2;
                        size_t const n{std::max(n1, n2)};
                        size_t const m{n / 2};
                        std::vector<uintmax_t> bits(m, 0);
                        std::copy(u_.bits().rbegin(),
                                  u_.bits().rbegin() + std::min(u_.bits().size(), m),
                                  bits.rbegin());
                        Integer const x0(bits);
                        bits = std::vector<uintmax_t>(m, 0);
                        std::copy(u_.bits().rbegin() + m,
                                  u_.bits().rbegin() + std::min(u_.bits().size(), 2 * m),
                                  bits.rbegin());
                        Integer const x1(bits);
                        bits = std::vector<uintmax_t>(m, 0);
                        std::copy(v_.bits().rbegin(),
                                  v_.bits().rbegin() + std::min(v_.bits().size(), m),
                                  bits.rbegin());
                        Integer const y0(bits);
                        bits = std::vector<uintmax_t>(m, 0);
                        std::copy(v_.bits().rbegin() + m,
                                  v_.bits().rbegin() + std::min(v_.bits().size(), 2 * m),
                                  bits.rbegin());
                        Integer const y1(bits);

                        auto const z0(IntegerMul<Integer, Integer>(x0, y0));
                        auto const z1(IntegerMul<Integer, Integer>(x1, y0)
                                      + IntegerMul<Integer, Integer>(x0, y1));
                        auto const z2(IntegerMul<Integer, Integer>(x1, y1));

                        //o = m * 8 * sizeof(uintmax_t)
                        //xy = z2 * 2^(2 * o) + z1 * 2^o + z0

                        Integer z(z0);
                        bits = std::vector<uintmax_t>(z1.bits().size() + m, 0);
                        std::copy(z1.bits().rbegin(), z1.bits().rend(), bits.rbegin() + m);
                        Integer const w1(bits);
                        z += w1;
                        bits = std::vector<uintmax_t>(z2.bits().size() + 2 * m, 0);
                        std::copy(z2.bits().rbegin(), z2.bits().rend(), bits.rbegin() + 2 * m);
                        Integer const w2(bits);
                        z += w2;

                        return z.bits();
                    }
                }
                else
                    return (u_ * -v_).bits();
            }
        }

        CONSTEXPR bool isNan() const
        {
            if (u_.isNan() || v_.isNan())
                return true;
            else
                return false;
        }

        CONSTEXPR bool isInfinity() const
        {
            return ((u_.isInfinity() && v_) || (u_ && v_.isInfinity()));
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerMul<E1, E2> operator*(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerMul<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
class IntegerDiv : public IntegerExpression<IntegerDiv<E1, E2> >
{
    public:
        CONSTEXPR IntegerDiv(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            if (u_.isPositive() && v_.isPositive())
                return true;
            else if (!u_.isPositive() && !v_.isPositive())
                return true;
            else
                return false;
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            if (!u_ || v_.isNan())
                return std::vector<uintmax_t>{};
            else if (v_.isInfinity())
                return std::vector<uintmax_t>{0};
            else
            {
                if (abs(Integer(u_)) < abs(Integer(v_)))
                    return std::vector<uintmax_t>{0};
                else if (u_.isPositive() && v_.isPositive())
                {
                    if (u_.template fits<uintmax_t>() && v_.template fits<uintmax_t>())
                        return std::vector<uintmax_t>{u_.template cast<uintmax_t>() / v_.template cast<uintmax_t>()};
                    else if (!(v_ & 1))
                    {
                        Integer r(v_);
                        Integer shift(0);

                        while (!(r & 1))
                        {
                            r >>= 1;
                            ++shift;
                        }

                        return ((u_ >> shift) / r).bits();
                    }
                    else
                        return computeQuotientBurnikelZiegler(u_, v_).bits();
                }
                else
                    return (u_ / -v_).bits();
            }
        }

        CONSTEXPR bool isNan() const
        {
            if (u_.isNan() || v_.isNan() || !v_)
                return true;
            else
                return false;
        }

        CONSTEXPR bool isInfinity() const
        {
            return (u_.isInfinity() && v_);
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerDiv<E1, E2> operator/(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerDiv<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E1, typename E2>
class IntegerMod : public IntegerExpression<IntegerMod<E1, E2> >
{
    public:
        CONSTEXPR IntegerMod(E1 const& u, E2 const& v) : u_(u), v_(v)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            if (u_.isPositive() && v_.isPositive())
                return true;
            else if (!u_.isPositive() && !v_.isPositive())
                return true;
            else
                return false;
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            if (!v_ || v_.isNan() || v_.isInfinity())
                return std::vector<uintmax_t>{};
            else
            {
                if ((u_.isPositive() && v_.isPositive()) ||
                    (!u_.isPositive() && !v_.isPositive()))
                {
                    if (v_ == 1)
                        return std::vector<uintmax_t>{0};
                    else if (v_ == 2)
                        return (u_ & 1).bits();
                    else if (abs(Integer(u_)).template fits<uintmax_t>() && abs(Integer(v_)).template fits<uintmax_t>())
                        return std::vector<uintmax_t>{abs(Integer(u_)).template cast<uintmax_t>()
                                                      % abs(Integer(v_)).template cast<uintmax_t>()};
                    else
                        return computeQrBurnikelZiegler(u_, v_).second.bits();
                }
                else
                    return computeQrBurnikelZiegler(u_, v_).second.bits();
            }
        }

        CONSTEXPR bool isNan() const
        {
            if (u_.isNan() || v_.isNan() || u_.isInfinity() || v_.isInfinity())
                return true;
            else
                return false;
        }

        CONSTEXPR bool isInfinity() const
        {
            return false;
        }

    private:
        E1 const u_;
        E2 const v_;
};

template <typename E1, typename E2>
CONSTEXPR inline IntegerMod<E1, E2> operator%(IntegerExpression<E1> const& u, IntegerExpression<E2> const& v)
{
    return IntegerMod<E1, E2>(*static_cast<const E1*>(&u), *static_cast<const E2*>(&v));
}

template <typename E>
class IntegerNeg : public IntegerExpression<IntegerNeg<E> >
{
    public:
        CONSTEXPR IntegerNeg(E const& u) : u_(u)
        {
        }

        CONSTEXPR bool isPositive() const
        {
            return !u_.isPositive();
        }

        CONSTEXPR std::vector<uintmax_t> bits() const
        {
            return u_.bits();
        }

        CONSTEXPR bool isNan() const
        {
            return u_.isNan();
        }

        CONSTEXPR bool isInfinity() const
        {
            return u_.isInfinity();
        }

    private:
        E const u_;
};

template <typename E>
CONSTEXPR inline IntegerNeg<E> operator-(IntegerExpression<E> const& u)
{
    return IntegerNeg<E>(*static_cast<const E*>(&u));
}

template <typename E>
class IntegerInv : public IntegerExpression<IntegerInv<E> >
{
public:
    CONSTEXPR IntegerInv(E const& u) : u_(u)
    {
    }

    CONSTEXPR bool isPositive() const
    {
        return !u_.isPositive();
    }

    CONSTEXPR std::vector<uintmax_t> bits() const
    {
        auto bits(u_.bits());

        auto threadFunc
            {
                [&bits] (size_t start, size_t end) -> void
                {
                    for (size_t i{start}; i < end; ++i)
                        bits[i] = ~bits[i];
                }
            };

        size_t const numThreads{std::thread::hardware_concurrency()};
        size_t const chunkSize{bits.size() / numThreads};
        std::vector<std::thread> threads;

        for (size_t i{0}; i < numThreads; ++i)
        {
            size_t const start{i * chunkSize};
            size_t const end{(i == numThreads - 1) ? bits.size() : (i + 1) * chunkSize};
            threads.emplace_back(threadFunc, start, end);
        }

        for (auto& t : threads)
            t.join();

        return bits;
    }

    CONSTEXPR bool isNan() const
    {
        return u_.isNan();
    }

    CONSTEXPR bool isInfinity() const
    {
        return u_.isInfinity();
    }

private:
    E const u_;
};

template <typename E>
CONSTEXPR inline IntegerInv<E> operator~(IntegerExpression<E> const& u)
{
    return IntegerInv<E>(*static_cast<const E*>(&u));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator+(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerAdd<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator+(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerAdd<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator-(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerSub<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator-(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerSub<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator*(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerMul<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator*(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerMul<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator/(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerDiv<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator/(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerDiv<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator%(IntegerExpression<E> const& lhs, T const& rhs)
{
    return IntegerMod<E, Integer >(*static_cast<const E *>(&lhs), Integer(rhs));
}

template <typename E, typename T, typename std::enable_if_t<std::is_standard_layout_v<T> && std::is_trivial_v<T> >* = nullptr>
CONSTEXPR inline decltype(auto) operator%(T const& lhs, IntegerExpression<E> const& rhs)
{
    return IntegerMod<Integer, E>(Integer(lhs), *static_cast<const E *>(&rhs));
}

inline std::pair<Integer, Integer> computeQrBurnikelZiegler(Integer const& dividend, Integer const& divisor)
{
    if (!divisor)
        return {Integer::nan(), Integer::nan()};
    else if (!dividend)
        return {Integer{0}, Integer{0}};
    else if (dividend.isNan())
        return {dividend, dividend};
    else if (divisor.isNan())
        return {divisor, divisor};
    else if (dividend.isInfinity() || divisor.isInfinity())
        return {Integer::nan(), Integer::nan()};
    else if (divisor.abs() > dividend.abs())
        return {Integer(0), dividend};
    else if (divisor.abs() == 1)
        return {Integer(divisor.sign() * dividend), Integer(0)};
    else if (dividend < 0 && divisor < 0)
    {
        auto qr{computeQrBurnikelZiegler(-dividend, -divisor)};

        if (qr.second)
        {
            ++qr.first;
            qr.second += divisor;
        }

        return qr;
    }
    else if (dividend > 0 && divisor < 0)
    {
        auto qr{computeQrBurnikelZiegler(dividend, -divisor)};

        qr.first = -qr.first;

        if (qr.second)
        {
            --qr.first;
            qr.second += divisor;
        }

        return qr;
    }
    else if (dividend < 0 && divisor > 0)
    {
        auto qr{computeQrBurnikelZiegler(-dividend, divisor)};

        qr.first = -qr.first;

        if (qr.second)
        {
            --qr.first;
            qr.second = divisor - qr.second;
        }

        return qr;
    }

    std::function<void(std::vector<Integer>&, Integer const&, Integer const&,
                       Integer const&, Integer const&)> inner1;

    inner1 =
    [&inner1] (std::vector<Integer>& a_digits, Integer const& x, Integer const& L,
               Integer const& R, Integer const& n) -> void
    {
        if (L + 1 == R)
        {
            a_digits[L] = x;
            return;
        }

        auto const mid((L + R) >> 1);
        auto const shift((mid - L) * n);
        auto const upper(x >> shift);
        auto const lower(x ^ (upper << shift));
        inner1(a_digits, lower, L, mid, n);
        inner1(a_digits, upper, mid, R, n);
    };

    auto _int2digits
    {
        [&inner1] (Integer const& a, Integer const& n) -> std::vector<Integer>
        {
            assert(a >= 0);

            if (!a)
                return std::vector<Integer>{Integer(0)};

            std::vector<Integer> a_digits(((a.number() + n - 1) / n), Integer(0));

            if (a)
                inner1(a_digits, a, Integer(0), Integer(a_digits.size()), n);

            return a_digits;
        }
    };

    std::function<Integer(std::vector<Integer> const&, Integer const&,
                          Integer const&, Integer const&)> inner2;

    inner2 =
    [&inner2] (std::vector<Integer> const& digits, Integer const& L,
               Integer const& R, Integer const& n) -> Integer
    {
        if (L + 1 == R)
            return digits[L];

        auto const mid((L + R) >> 1);
        auto const shift((mid - L) * n);

        return (inner2(digits, mid, R, n) << shift) + inner2(digits, L, mid, n);
    };

    auto _digits2int
    {
        [&inner2] (std::vector<Integer> const& digits, Integer const& n) -> Integer
        {
            if (digits.empty())
                return Integer(0);

            return inner2(digits, Integer(0), Integer(digits.size()), n);
        }
    };

    std::function<std::pair<Integer, Integer>(Integer, Integer, Integer)> _div2n1n;
    std::function<std::pair<Integer, Integer>(Integer const&, Integer const&,
                                              Integer const&, Integer const&,
                                              Integer const&, Integer const&)> _div3n2n;

    _div2n1n =
    [&_div3n2n] (Integer a, Integer b, Integer n) -> std::pair<Integer, Integer>
    {
        if (a.fits<uintmax_t>() && b.fits<uintmax_t>())
            return {Integer(a.cast<uintmax_t>() / b.cast<uintmax_t>()),
                    Integer(a.cast<uintmax_t>() % b.cast<uintmax_t>())};

        auto pad(n & 1);

        if (pad)
        {
            a <<= 1;
            b <<= 1;
            ++n;
        }

        auto const half_n(n >> 1);
        auto const mask((Integer(1) << half_n) - 1);
        auto const b1(b >> half_n);
        auto const b2(b & mask);
        auto[q1, r] = _div3n2n(a >> n, (a >> half_n) & mask, b, b1, b2, half_n);
        auto[q2, r2] = _div3n2n(r, a & mask, b, b1, b2, half_n);
        r = r2;

        if (pad)
            r >>= 1;

        return {(q1 << half_n) | q2, r};
    };

    _div3n2n =
    [&_div2n1n] (Integer const& a12, Integer const& a3,
                 Integer const& b, Integer const& b1,
                 Integer const& b2, Integer const& n) -> std::pair<Integer, Integer>
    {
        Integer q, r;

        if (a12 >> n == b1)
        {
            q = (Integer(1) << n) - 1;
            r = a12 - (b1 << n) + b1;
        }
        else
        {
            auto const p{_div2n1n(a12, b1, n)};
            q = p.first;
            r = p.second;
        }

        r = ((r << n) | a3) - q * b2;

        while (r < 0)
        {
            --q;
            r += b;
        }

        return {q, r};
    };

    auto const n{divisor.number()};
    auto const a_digits(_int2digits(dividend, n));

    Integer r(0);
    Integer q(0);
    std::vector<Integer> q_digits;

    for (auto it{a_digits.rbegin()}; it != a_digits.rend(); ++it)
    {
        auto[q_digit, r_] = _div2n1n((r << n) + *it, divisor, n);
        r = r_;
        q_digits.emplace_back(q_digit);
    }

    std::reverse(q_digits.begin(), q_digits.end());

    q = _digits2int(q_digits, n);

    return {q, r};
}

inline Integer computeQuotientBurnikelZiegler(Integer const& dividend, Integer const& divisor)
{
    if (!divisor)
        return Integer::nan();
    else if (!dividend)
        return Integer{0};
    else if (dividend.isNan())
        return dividend;
    else if (divisor.isNan())
        return divisor;
    else if (dividend.isInfinity() || divisor.isInfinity())
        return Integer::nan();
    else if (divisor.abs() > dividend.abs())
        return Integer{0};

    auto qr{computeQrBurnikelZiegler(dividend, divisor)};

    if (dividend < 0 && divisor < 0)
    {
        if (qr.second)
            --qr.first;
    }
    else if (dividend > 0 && divisor < 0)
    {
        if (qr.second)
            ++qr.first;
    }
    else if (dividend < 0 && divisor > 0)
    {
        if (qr.second)
            ++qr.first;
    }

    return qr.first;
}

template <typename S>
CONSTEXPR inline std::pair<Integer, Integer> computeQrBurnikelZiegler(Integer const& dividend, S const& divisor)
{
    return computeQrBurnikelZiegler(dividend, Integer(divisor));
}

template <typename S>
CONSTEXPR inline std::pair<Integer, Integer> computeQrBurnikelZiegler(S const& dividend, Integer const& divisor)
{
    return computeQrBurnikelZiegler(Integer(dividend), divisor);
}

#endif // INTEGER_H
