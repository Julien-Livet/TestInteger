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

        CONSTEXPR bool isPositive() const noexcept
        {
            return static_cast<E const&>(*this).isPositive();
        }

        CONSTEXPR bool isNegative() const noexcept
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

        template <size_t N>
        CONSTEXPR Integer(std::bitset<N> const& bits, bool isPositive = true) : isPositive_{isPositive}
        {
            setBits(0, bits);
        }

        template <class InputIt>
        CONSTEXPR Integer(InputIt begin, InputIt end, bool isPositive = true) : isPositive_{isPositive}, bits_{begin, end}
        {
            if (bits_.empty())
                bits_.emplace_back(0);

            adjust();
        }

        Integer(std::vector<uintmax_t> const& bits, bool isPositive = true);
        Integer(std::initializer_list<uintmax_t> const& bits, bool isPositive = true);

#ifdef WITH_GMP
        Integer(mpz_class const& n);
#endif

        Integer(char const* n, size_t base = 0);
        Integer(std::string n, size_t base = 0);
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
        std::string toString(size_t base = 10) const;
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
        size_t precision() const noexcept;
        void setPrecision(size_t precision);

        template <typename URNG>
        CONSTEXPR void setRandom()
        {
            isNan_ = false;
            isInfinity_ = false;

            URNG g;

            isPositive_ = g() % 2;

            auto threadFunc
            {
                [this] (size_t start, size_t end) -> void
                {
                    URNG g;

                    for (size_t i{start}; i < end; ++i)
                    {
                        auto const n{g()};

                        if (sizeof(uintmax_t) <= sizeof(n))
                            this->bits_[i] = static_cast<uintmax_t>(n);
                        else
                        {
                            this->bits_[i] = 0;

                            auto const jMax{sizeof(uintmax_t) / sizeof(n)};

                            for (size_t j{0}; j < jMax; ++j)
                            {
                                this->bits_[i] <<= sizeof(n) * 8;
                                this->bits_[i] |= g();
                            }
                        }
                    }
                }
            };

            size_t const numThreads{std::thread::hardware_concurrency()};
            size_t const chunkSize{bits_.size() / numThreads};
            std::vector<std::thread> threads;

            for (size_t i{0}; i < numThreads; ++i)
            {
                size_t const start{i * chunkSize};
                size_t const end{(i == numThreads - 1) ? bits_.size() : (i + 1) * chunkSize};
                threads.emplace_back(threadFunc, start, end);
            }

            for (auto& t : threads)
                t.join();
        }

        int isPrime(size_t reps = 50) const;
        CONSTEXPR void setPositive();
        CONSTEXPR void setNegative();

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

        bool bit(size_t n) const noexcept;
        void setBit(size_t n, bool bit);
        uintmax_t bits(size_t n) const noexcept;
        void setBits(size_t n, uintmax_t const& bits);

        template <size_t N>
        CONSTEXPR void setBits(size_t n, std::bitset<N> const& bits)
        {
            for (size_t i{0}; i < bits.size(); ++i)
                setBit(n + i, bits[i]);

            if (autoAdjust_)
                adjust();
        }

        size_t count() const noexcept;
        Integer number() const noexcept;
        bool isEven() const noexcept;
        bool isOdd() const noexcept;

        template <typename S>
        CONSTEXPR bool fits() const
        {
            return (*this == this->cast<S>());
        }

        Integer sign() const;
        CONSTEXPR void setSign(Integer const& other) noexcept;
        Integer previousPrime() const;
        Integer nextPrime() const;
        size_t size() const noexcept;
        void adjust();
        bool isCoprime(Integer const& other) const noexcept;
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

inline Integer const& min(Integer const& a, Integer const& b)
{
    return a < b ? a : b;
}

inline Integer const& max(Integer const& a, Integer const& b)
{
    return a > b ? a : b;
}

inline Integer number(Integer const& n) noexcept
{
    return n.number();
}

#ifdef WITH_GMP
template <>
inline mpz_class Integer::cast<mpz_class>() const
{
    auto s{toString(2)};

    s.replace(s.find("0b"), 2, "");

    return mpz_class{s, 2};
}
#endif

template <>
inline std::string Integer::cast<std::string>() const
{
    return toString();
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

inline std::ostream& operator<<(std::ostream& os, Integer const& n)
{
    return os << n.toString();
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

inline Integer gcd(Integer const& a, Integer const& b)
{
    if (a.isNan() || b.isNan() || a.isInfinity() || b.isInfinity())
        return Integer::nan();

    if (a < 0)
        return gcd(a.abs(), b);

    if (b < 0)
        return gcd(a, b.abs());

    if (a < b)
        return gcd(b, a);

    if (!a)
        return b;

    if (!b)
        return a;

    if (a.isEven() && b.isEven())
        return 2 * gcd(a >> 1, b >> 1);
    else if (a.isOdd() && b.isEven())
        return gcd(a, Integer(b >> 1));
    else if (a.isEven() && b.isOdd())
        return gcd(Integer(a >> 1), b);
    else //if (a.isOdd() && b.isOdd())
        return gcd(Integer((a - b) >> 1), b);
}

template <typename S>
CONSTEXPR inline Integer gcd(Integer const& a, S const& b)
{
    return gcd(a, Integer(b));
}

template <typename S>
CONSTEXPR inline Integer gcd(S const& a, Integer const& b)
{
    return gcd(Integer(a), b);
}

inline Integer lcm(Integer const& a, Integer const& b)
{
    return abs(Integer(a * b)) / gcd(a, b);
}

template <typename S>
CONSTEXPR inline Integer lcm(Integer const& a, S const& b)
{
    return lcm(a, Integer(b));
}

template <typename S>
CONSTEXPR inline Integer lcm(S const& a, Integer const& b)
{
    return lcm(Integer(a), b);
}

inline Integer gcdExtended(Integer a, Integer b, Integer& u, Integer& v)
{
    if (a.isNan() || b.isNan() || a.isInfinity() || b.isInfinity())
    {
        u.setNan();
        v.setNan();

        return Integer::nan();
    }

    if (!a && !b)
        return Integer(0);

    Integer r1(a), u1(1), v1(0);
    Integer r2(b), u2(0), v2(1);
    Integer q, r_temp, u_temp, v_temp;

    while (r2 != 0)
    {
        q = r1 / r2;
        r_temp = r1 - q * r2;
        u_temp = u1 - q * u2;
        v_temp = v1 - q * v2;

        r1 = r2;
        u1 = u2;
        v1 = v2;
        r2 = r_temp;
        u2 = u_temp;
        v2 = v_temp;
    }

    u = u1;
    v = v1;

    return r1;
}

template <typename S1, typename S2>
CONSTEXPR inline Integer gcdExtended(S1 const& a, S2 const& b, Integer& u, Integer& v)
{
    return gcdExtended(Integer(a), Integer(b), u, v);
}

inline Integer pow(Integer base, Integer exp)
{
    assert(exp >= 0);

    if (base.isInfinity() || base.isNan())
        return base;
    else if (exp.isNan() || exp.isInfinity())
        return exp;
    else if (base < 0)
    {
        auto n(pow(base.abs(), exp));

        if (exp & 1)
            n = -n;

        return n;
    }
    else if (base == 2)
        return Integer(1) << exp;

    Integer result(1);

    for (;;)
    {
        if (exp & 1)
            result *= base;

        exp >>= 1;

        if (!exp)
            break;

        base *= base;
    }

    return result;
}

template <typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline Integer pow(Integer const& base, S const& exp)
{
    return pow(base, Integer(exp));
}

template <typename S, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr>
CONSTEXPR inline Integer pow(S const& base, Integer const& exp)
{
    return pow(Integer(base), exp);
}

inline Integer factorial(Integer const& n)
{
    if (n.isNan() || n.isInfinity())
        return n;

    assert(n >= 0);

    if (n == 0)
        return Integer(1);

    return n * factorial(n - 1);
}

inline Integer doubleFactorial(Integer const& n)
{
    return factorial(factorial(n));
}

inline Integer multiFactorial(Integer const& n, Integer const& m)
{
    return pow(factorial(n), m);
}

template <typename S>
CONSTEXPR inline Integer multiFactorial(Integer const& n, S const& m)
{
    return multiFactorial(n, Integer(m));
}

template <typename S>
CONSTEXPR inline Integer multiFactorial(S const& n, Integer const& m)
{
    return multiFactorial(Integer(n), m);
}

inline Integer powm(Integer base, Integer exp, Integer const& mod)
{
    assert(exp >= 0);

    Integer result(1);
    Integer base_mod(base % mod);

    while (exp > 0)
    {
        if ((exp & 1) == 1)
        {
            result *= base_mod;
            result %= mod;
        }

        base_mod *= base_mod;
        base_mod %= mod;

        exp >>= 1;
    }

    return result;
}

template <typename S, typename U, typename std::enable_if_t<std::is_standard_layout_v<S> && std::is_trivial_v<S> >* = nullptr, typename std::enable_if_t<std::is_standard_layout_v<U> && std::is_trivial_v<U> >* = nullptr>
CONSTEXPR inline Integer powm(Integer const& base, S const& exp, U const& mod)
{
    return powm(base, Integer(exp), Integer(mod));
}

inline std::pair<Integer, Integer> computeQr(Integer const& dividend, Integer const& divisor)
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
        return {Integer{0}, dividend};
    else if (dividend < 0 && divisor < 0)
    {
        auto qr{computeQr(-dividend, -divisor)};

        if (qr.second)
        {
            ++qr.first;
            qr.second += divisor;
        }

        return qr;
    }
    else if (dividend > 0 && divisor < 0)
    {
        auto qr{computeQr(dividend, -divisor)};

        qr.first = -qr.first;

        if (qr.second)
        {
            qr.first -= 1;
            qr.second += divisor;
        }

        return qr;
    }
    else if (dividend < 0 && divisor > 0)
    {
        auto qr{computeQr(-dividend, divisor)};

        qr.first = -qr.first;

        if (qr.second)
        {
            qr.first -= 1;
            qr.second = divisor - qr.second;
        }

        return qr;
    }

    Integer start(1);
    auto end(dividend);

    while (start <= end)
    {
        Integer mid(end + start);
        mid >>= 1;

        Integer n(dividend - divisor * mid);

        if (n > divisor)
            start = mid + 1;
        else if (n < 0)
            end = mid - 1;
        else
        {
            if (n == divisor)
            {
                ++mid;
                n = 0;
            }

            return {Integer(mid), n};
        }
    }

    return {Integer(0), dividend};
}

template <typename S>
CONSTEXPR inline std::pair<Integer, Integer> computeQr(Integer const& dividend, S const& divisor)
{
    return computeQr(dividend, Integer(divisor));
}

template <typename S>
CONSTEXPR inline std::pair<Integer, Integer> computeQr(S const& dividend, Integer const& divisor)
{
    return computeQr(Integer(dividend), divisor);
}

inline Integer computeQuotient(Integer const& dividend, Integer const& divisor)
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
    else if (dividend < 0 && divisor < 0)
        return computeQuotient(-dividend, -divisor);
    else if (dividend > 0 && divisor < 0)
        return -computeQuotient(dividend, -divisor);
    else if (dividend < 0 && divisor > 0)
        return -computeQuotient(-dividend, divisor);

    Integer start(1);
    auto end(dividend);

    while (start <= end)
    {
        Integer mid(end + start);
        mid >>= 1;

        Integer n(dividend - divisor * mid);

        if (n > divisor)
            start = mid + 1;
        else if (n < 0)
            end = mid - 1;
        else
        {
            if (n == divisor)
            {
                ++mid;
                n = 0;
            }

            return mid;
        }
    }

    return Integer(0);
}

template <typename S>
CONSTEXPR inline Integer computeQuotient(Integer const& dividend, S const& divisor)
{
    return computeQuotient(dividend, Integer(divisor));
}

template <typename S>
CONSTEXPR inline Integer computeQuotient(S const& dividend, Integer const& divisor)
{
    return computeQuotient(Integer(dividend), divisor);
}

template <typename S>
CONSTEXPR inline Integer fibonacci(Integer n)
{
    assert(n >= 0);

    if (!n)
        return 0;
    else if (n == 1)
        return 1;

    n -= 1;

    Integer fn_2(0);
    Integer fn_1(1);
    Integer fn;

    while (n)
    {
        fn = fn_1 + fn_2;
        fn_2 = fn_1;
        fn_1 = fn;
        --n;
    }

    return fn;
}

inline Integer primorial(Integer n)
{
    Integer result(1);
    Integer number(2);

    while (number <= n)
    {
        result *= number;
        number = number.nextPrime();
    }

    return result;
}

inline int legendre(Integer const& a, Integer const& p)
{
    assert(p.isPrime());

    if (!(a % p))
        return 0;
    else
    {
        bool isResidue{false};

        if (p == 2)
            isResidue = true;
        else
            isResidue = (powm(a, (p - 1) / 2, p) == 1);

        if (isResidue)
            return 1;
        else
            return -1;
    }
}

template <typename S>
CONSTEXPR inline Integer legendre(Integer const& a, S const& p)
{
    return legendre(a, Integer(p));
}

template <typename S>
CONSTEXPR inline Integer legendre(S const& a, Integer const& p)
{
    return legendre(Integer(a), p);
}

inline int jacobi(Integer const& a, Integer n)
{
    assert(n > 0 && n.isOdd());

    int result(1);
    Integer prime(2);

    while (n != 1)
    {
        if (!(n % prime))
        {
            n /= prime;
            result *= legendre(a, prime);
        }
        else
            prime = prime.nextPrime();
    }

    return result;
}

template <typename S>
CONSTEXPR inline int jacobi(Integer const& a, S const& n)
{
    return jacobi(a, Integer(n));
}

template <typename S>
CONSTEXPR inline int jacobi(S const& a, Integer const& n)
{
    return jacobi(Integer(a), n);
}

CONSTEXPR inline int kronecker(Integer const& a, Integer const& b)
{
    if (a == b)
        return 1;
    else
        return 0;
}

template <typename S>
CONSTEXPR inline int kronecker(Integer const& a, S const& b)
{
    return kronecker(a, Integer(b));
}

template <typename S>
CONSTEXPR inline int kronecker(S const& a, Integer const& b)
{
    return kronecker(Integer(a), b);
}

inline Integer binomial(Integer const& n, Integer const& k)
{
    assert(n >= 0 && k >= 0);

    return factorial(n) / (factorial(k) * factorial(n - k));
}

template <typename S>
CONSTEXPR inline Integer binomial(Integer const& n, S const& k)
{
    return binomial(n, Integer(k));
}

template <typename S>
CONSTEXPR inline Integer binomial(S const& n, Integer const& k)
{
    return binomial(Integer(n), k);
}

inline Integer sqrt(Integer const& n)
{
    if (n < 0)
        return Integer::nan();
    else if (!n || n == 1 || n.isNan() || n.isInfinity())
        return n;

    Integer lo(1), hi(n);
    Integer res(1);

    while (lo <= hi)
    {
        Integer mid(lo + hi);
        mid >>= 1;

        if (mid * mid <= n)
        {
            res = mid;
            lo = mid + 1;
        }
        else
            hi = mid - 1;
    }

    return res;
}

inline Integer root(Integer const& x, Integer const& n)
{
    assert(n > 0);

    if (x < 0)
        return Integer::nan();
    else if (!x || x == 1 || x.isNan() || x.isInfinity())
        return x;
    else if (n == 1)
        return x;
    else if (n.isNan() || n.isInfinity())
        return Integer::nan();

    Integer lo(1), hi(x);
    Integer res(1);

    while (lo <= hi)
    {
        Integer mid(lo + hi);
        mid >>= 1;

        if (pow(mid, n) <= x)
        {
            res = mid;
            lo = mid + 1;
        }
        else
            hi = mid - 1;
    }

    return res;
}

template <typename S>
CONSTEXPR inline Integer root(S const& x, Integer const& n)
{
    return root(Integer(x), n);
}

template <typename S>
CONSTEXPR inline Integer root(Integer const& x, S const& n)
{
    return root(x, Integer(n));
}

inline Integer computeQuotientBinary(Integer dividend, Integer const& divisor)
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
    else if (dividend < 0 && divisor > 0)
        return -computeQuotientBinary(-dividend, divisor);
    else if (dividend > 0 && divisor < 0)
        return -computeQuotientBinary(dividend, -divisor);
    else if (dividend < 0 && divisor < 0)
        return computeQuotientBinary(-dividend, -divisor);

    Integer quotient(0);
    auto tempDivisor(divisor);
    Integer bit(1);

    while (dividend >= (tempDivisor << 1))
    {
        tempDivisor <<= 1;
        bit <<= 1;
    }

    while (bit >= 1)
    {
        if (dividend >= tempDivisor)
        {
            dividend -= tempDivisor;
            quotient += bit;
        }

        tempDivisor >>= 1;
        bit >>= 1;
    }

    return quotient;
}

template <typename S>
CONSTEXPR inline Integer computeQuotientBinary(Integer const& dividend, S const& divisor)
{
    return computeQuotientBinary(dividend, Integer(divisor));
}

template <typename S>
CONSTEXPR inline Integer computeQuotientBinary(S const& dividend, Integer const& divisor)
{
    return computeQuotientBinary(Integer(dividend), divisor);
}

inline std::pair<Integer, Integer> computeQrBinary(Integer const& dividend, Integer const& divisor)
{
    std::pair<Integer, Integer> qr{computeQuotientBinary(dividend, divisor), Integer(0)};
    qr.second = dividend.abs() - qr.first.abs() * divisor.abs();

    if (dividend < 0 && divisor < 0)
    {
        if (qr.second)
        {
            ++qr.first;
            qr.second += divisor;
        }
    }
    else if (dividend > 0 && divisor < 0)
    {
        if (qr.second)
        {
            qr.first -= 1;
            qr.second += divisor;
        }
    }
    else if (dividend < 0 && divisor > 0)
    {
        if (qr.second)
        {
            qr.first -= 1;
            qr.second = divisor - qr.second;
        }
    }

    return qr;
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

inline Integer operator""_z(char const* str)
{
    return Integer(str);
}

#endif // INTEGER_H
