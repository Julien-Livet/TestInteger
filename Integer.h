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
#include <iomanip>
#include <random>
#include <sstream>
#include <vector>

#if __cplusplus >= 201703L
#define CONSTEXPR constexpr
#else
#define CONSTEXPR
#endif

#if __cplusplus >= 202002L
#include <format>
#endif

#ifdef USING_GMP
#include <gmpxx.h>
#endif

using longest_type = uintmax_t;

template <typename T, typename Enable = void>
class Integer;

template <typename T>
class Integer<T, typename std::enable_if<std::is_unsigned<T>::value>::type>
{
    public:
        CONSTEXPR Integer() : bits_{T{0}}
        {
        }

        template <typename S>
        CONSTEXPR explicit Integer(S n) : isPositive_{n >= 0}
        {
            bits_.reserve(std::max(longest_type{1}, longest_type{sizeof(S) / sizeof(T)}));

            if (n < 0)
                n = -n;

            if (sizeof(T) == sizeof(S))
                bits_.emplace_back(n);
            else
            {
                auto const shift{longest_type{1} << std::min(sizeof(T), sizeof(S)) * 8};

                for (size_t i{0}; i < bits_.capacity(); ++i)
                {
                    bits_.emplace_back(n % shift);
                    n /= shift;
                }

                std::reverse(bits_.begin(), bits_.end());
            }
        }

        CONSTEXPR explicit Integer(std::vector<T> const& bits, bool isPositive = true) : isPositive_{isPositive}, bits_{bits}
        {
            if (bits_.empty())
                bits_.emplace_back(T{0});

            adjust();
        }

        template <size_t N>
        CONSTEXPR explicit Integer(std::bitset<N> const& bits, bool isPositive = true) : isPositive_{isPositive}
        {
            setBits(0, bits);
        }

        CONSTEXPR explicit Integer(std::initializer_list<T> const& bits, bool isPositive = true) : isPositive_{isPositive}, bits_{bits}
        {
            if (bits_.empty())
                bits_.emplace_back(T{0});

            adjust();
        }

        template <class InputIt>
        CONSTEXPR explicit Integer(InputIt begin, InputIt end, bool isPositive = true) : isPositive_{isPositive}, bits_{begin, end}
        {
            if (bits_.empty())
                bits_.emplace_back(T{0});

            adjust();
        }

#ifdef USING_GMP
        CONSTEXPR explicit Integer(mpz_class const& n) : Integer(n.get_str(2), 2)
        {
        }
#endif

        CONSTEXPR explicit Integer(char const* n, size_t base = 0) : Integer(std::string{n}, base)
        {
        }

        CONSTEXPR explicit Integer(std::string n, size_t base = 0)
        {
            n.erase(std::remove_if(n.begin(), n.end(), isspace), n.end());
            n.erase(std::remove(n.begin(), n.end(), '\''), n.end());

            auto it{n.begin()};

            if (*it == '-')
            {
                isPositive_ = false;
                ++it;
            }

            if (!base)
            {
                auto s{std::string{it, n.end()}.substr(0, 2)};
                std::transform(s.begin(), s.end(), s.begin(),
                               [] (unsigned char c) { return std::tolower(c); });

                if (s[0] == 'b' || s == "0b")
                    base = 2;
                else if (s[0] == 'o' || s == "0o")
                    base = 8;
                else if (s[0] == 'x' || s == "0x")
                    base = 16;
                else
                    base = 10;
            }

            assert(2 <= base && base <= 62);

            std::string str{it, n.end()};
            std::transform(str.begin(), str.end(), str.begin(),
                           [] (unsigned char c) { return std::tolower(c); });

            if (str == "nan")
                setNan();
            else if (str == "inf")
                setInfinity();
            else
            {
                auto const isPositive{isPositive_};

                if (base == 2)
                {
                    if (*it == 'b' || *it == 'B')
                        ++it;
                    else if (std::string(it, n.end()).substr(0, 2) == "0b"
                             || std::string(it, n.end()).substr(0, 2) == "0B")
                        it += 2;

                    while (it != n.end())
                    {
                        if (*it == '1')
                        {
                            *this |= 1;

                            if (it != n.end() - 1)
                                *this <<= 1;
                        }
                        else if (*it == '0')
                        {
                            if (it != n.end() - 1)
                                *this <<= 1;
                        }

                        ++it;
                    }
                }
                else if (base == 8)
                {
                    if (*it == 'o' || *it == 'O')
                        ++it;
                    else if (std::string(it, n.end()).substr(0, 2) == "0o"
                             || std::string(it, n.end()).substr(0, 2) == "0O")
                        it += 2;

                    auto otherIt{n.rbegin()};
                    Integer p(1);

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '7')
                        {
                            *this += (*otherIt - '0') * p;
                            p *= base;
                        }

                        ++otherIt;
                    }
                }
                else if (base <= 10)
                {
                    auto otherIt{n.rbegin()};
                    Integer p(1);

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= static_cast<char>('0' + base))
                        {
                            *this += (*otherIt - '0') * p;
                            p *= base;
                        }

                        ++otherIt;
                    }
                }
                else if (base < 16)
                {
                    auto otherIt{n.rbegin()};
                    Integer p(1);

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '9')
                        {
                            *this += (*otherIt - '0') * p;
                            p *= base;
                        }
                        else if ('a' <= std::tolower(*otherIt) && std::tolower(*otherIt) <= static_cast<char>('a' + base - 10))
                        {
                            *this += (*otherIt - 'a' + 10) * p;
                            p *= base;
                        }

                        ++otherIt;
                    }
                }
                else if (base == 16)
                {
                    if (*it == 'x' || *it == 'X')
                        ++it;
                    else if (std::string(it, n.end()).substr(0, 2) == "0x"
                             || std::string(it, n.end()).substr(0, 2) == "0X")
                        it += 2;

                    auto otherIt{n.rbegin()};
                    Integer p(1);

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '9')
                        {
                            *this += (*otherIt - '0') * p;
                            p *= base;
                        }
                        else if ('a' <= std::tolower(*otherIt) && std::tolower(*otherIt) <= 'f')
                        {
                            *this += (*otherIt - 'a' + 10) * p;
                            p *= base;
                        }

                        ++otherIt;
                    }
                }
                else// if (base <= 62)
                {
                    auto otherIt{n.rbegin()};
                    Integer p(1);

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '9')
                        {
                            *this += (*otherIt - '0') * p;
                            p *= base;
                        }
                        else if ('a' <= *otherIt && *otherIt <= 'z')
                        {
                            *this += (*otherIt - 'a' + 10) * p;
                            p *= base;
                        }
                        else if ('A' <= *otherIt && *otherIt <= 'Z')
                        {
                            *this += (*otherIt - 'A' + 36) * p;
                            p *= base;
                        }

                        ++otherIt;
                    }
                }

                isPositive_ = isPositive;
            }

            adjust();
        }

        template <typename S>
        CONSTEXPR explicit Integer(Integer<S> const& other) : Integer(other.toString(2), 2)
        {
        }

        CONSTEXPR bool isPositive() const noexcept
        {
            return isPositive_;
        }

        CONSTEXPR bool isNegative() const noexcept
        {
            return !isPositive_;
        }

        CONSTEXPR auto const& bits() const noexcept
        {
            return bits_;
        }

        CONSTEXPR void invert() noexcept
        {
            for (size_t i{0}; i < bits_.size(); ++i)
                bits_[i] = ~bits_[i];

            if (autoAdjust_)
                adjust();
        }

        CONSTEXPR Integer& operator*=(Integer const& other)
        {
            auto const lhs(*this);
            auto const rhs(other);

            if (other.isNegative())
            {
                *this = -*this;

                return *this *= -other;
            }

            if (isNan() || other.isNan())
                setNan();
            else if (!*this || !other)
            {
                *this = 0;
            }
            else if (isInfinity() || other.isInfinity())
            {
                setInfinity();

                if (other.isNegative())
                    isPositive_ = !isPositive_;
            }
            else if (!*this || !other)
                *this = 0;
            else
            {
                if (isPositive_ && other.isPositive_)
                {
                    if (this->template fits<longest_type>() && other.template fits<longest_type>())
                    {
                        auto const a{this->template cast<longest_type>()};
                        auto const b{other.template cast<longest_type>()};
                        auto const ab{a * b};

                        if (ab / b == a)
                            *this = ab;
                        else
                        {
                            auto number{[] (T n) -> size_t
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
                            size_t n{std::max(number(bits_.back()), number(other.bits_.back()))};
                            if (n % 2)
                                ++n;
                            size_t const m{n / 2};
                            Integer x0, x1, y0, y1;
                            x0.bits_ = std::vector<T>(1, T{0});
                            x0.bits_.back() = (static_cast<T>(~T{0}) >> (sizeof(T) * 8 - m)) & bits_.back();
                            x1.bits_ = std::vector<T>(1, T{0});
                            x1.bits_.back() = (((static_cast<T>(~T{0}) >> (sizeof(T) * 8 - m)) << m) & bits_.back()) >> m;
                            y0.bits_ = std::vector<T>(1, T{0});
                            y0.bits_.back() = (static_cast<T>(~T{0}) >> (sizeof(T) * 8 - m)) & other.bits_.back();
                            y1.bits_ = std::vector<T>(1, T{0});
                            y1.bits_.back() = (((static_cast<T>(~T{0}) >> (sizeof(T) * 8 - m)) << m) & other.bits_.back()) >> m;

                            assert(*this == ((x1 << m) | x0));
                            assert(other == ((y1 << m) | y0));

                            auto const z0(x0 * y0);
                            auto const z1(x1 * y0 + x0 * y1);
                            auto const z2(x1 * y1);

#ifdef USING_GMP
                            assert(z0 == mpz_class{x0.template cast<mpz_class>() * y0.template cast<mpz_class>()});
                            assert(z1 == mpz_class{x1.template cast<mpz_class>() * y0.template cast<mpz_class>() + x0.template cast<mpz_class>() * y1.template cast<mpz_class>()});
                            assert(z2 == mpz_class{x1.template cast<mpz_class>() * y1.template cast<mpz_class>()});

                            mpz_class _2_2m{2};
                            mpz_pow_ui(_2_2m.get_mpz_t(), _2_2m.get_mpz_t(), 2 * m);
                            mpz_class _2_m{2};
                            mpz_pow_ui(_2_m.get_mpz_t(), _2_m.get_mpz_t(), m);
                            mpz_class const n1_{lhs.template cast<mpz_class>() * rhs.template cast<mpz_class>()};
                            mpz_class const n2_{z2.template cast<mpz_class>() * _2_2m
                                                + z1.template cast<mpz_class>() * _2_m
                                                + z0.template cast<mpz_class>()};
                            assert(n1_ == n2_);
#endif

                            //xy = z2 * 2^(2 * m) + z1 * 2^m + z0
                            *this = z0 + (z1 << m) + (z2 << 2 * m);
                        }
                    }
                    else if (!(rhs & 1))
                    {
                        auto r(rhs);
                        Integer shift(0);

                        while (!(r & 1))
                        {
                            r >>= 1;
                            ++shift;
                        }

                        *this <<= shift;
                        *this *= r;
                    }
                    else
                    {
                        //Karatsuba algorithm
                        //x = x1 * 2^m + x0
                        //y = y1 * 2^m + y0
                        size_t n1{number() / (sizeof(T) * 8)};
                        if (number() % (sizeof(T) * 8))
                            ++n1;
                        if (n1 % 2)
                            ++n1;
                        size_t n2{other.number() / (sizeof(T) * 8)};
                        if (other.number() % (sizeof(T) * 8))
                            ++n2;
                        if (n2 % 2)
                            ++n2;
                        size_t const n{std::max(n1, n2)};
                        size_t const m{n / 2};
                        Integer x0, x1, y0, y1;
                        x0.bits_ = std::vector<T>(m, T{0});
                        std::copy(bits_.rbegin(),
                                  bits_.rbegin() + std::min(bits_.size(), m),
                                  x0.bits_.rbegin());
                        x1.bits_ = std::vector<T>(m, T{0});
                        std::copy(bits_.rbegin() + m,
                                  bits_.rbegin() + std::min(bits_.size(), 2 * m),
                                  x1.bits_.rbegin());
                        y0.bits_ = std::vector<T>(m, T{0});
                        std::copy(other.bits_.rbegin(),
                                  other.bits_.rbegin() + std::min(other.bits_.size(), m),
                                  y0.bits_.rbegin());
                        y1.bits_ = std::vector<T>(m, T{0});
                        std::copy(other.bits_.rbegin() + m,
                                  other.bits_.rbegin() + std::min(other.bits_.size(), 2 * m),
                                  y1.bits_.rbegin());

                        assert(*this == ((x1 << (m * sizeof(T) * 8)) | x0));
                        assert(other == ((y1 << (m * sizeof(T) * 8)) | y0));

                        auto const z0(x0 * y0);
                        auto const z1(x1 * y0 + x0 * y1);
                        auto const z2(x1 * y1);

#ifdef USING_GMP
                        assert(z0 == mpz_class{x0.template cast<mpz_class>() * y0.template cast<mpz_class>()});
                        assert(z1 == mpz_class{x1.template cast<mpz_class>() * y0.template cast<mpz_class>() + x0.template cast<mpz_class>() * y1.template cast<mpz_class>()});
                        assert(z2 == mpz_class{x1.template cast<mpz_class>() * y1.template cast<mpz_class>()});
#endif

                        //o = m * 8 * sizeof(T)
                        //xy = z2 * 2^(2 * o) + z1 * 2^o + z0

#ifdef USING_GMP
                        size_t const o{m * 8 * sizeof(T)};
                        mpz_class _2_2o{2};
                        mpz_pow_ui(_2_2o.get_mpz_t(), _2_2o.get_mpz_t(), 2 * o);
                        mpz_class _2_o{2};
                        mpz_pow_ui(_2_o.get_mpz_t(), _2_o.get_mpz_t(), o);
                        mpz_class const n1_{lhs.template cast<mpz_class>() * rhs.template cast<mpz_class>()};
                        mpz_class const n2_{z2.template cast<mpz_class>() * _2_2o
                                            + z1.template cast<mpz_class>() * _2_o
                                            + z0.template cast<mpz_class>()};
                        assert(n1_ == n2_);
#endif

                        *this = z0;
                        Integer w1, w2;
                        w1.bits_ = std::vector<T>(z1.bits_.size() + m, T{0});
                        std::copy(z1.bits_.rbegin(), z1.bits_.rend(), w1.bits_.rbegin() + m);
                        *this += w1;
                        w2.bits_ = std::vector<T>(z2.bits_.size() + 2 * m, T{0});
                        std::copy(z2.bits_.rbegin(), z2.bits_.rend(), w2.bits_.rbegin() + 2 * m);
                        *this += w2;
                    }
                }
                else
                {
                    *this *= -other;
                    *this = -*this;
                }
            }

            if (autoAdjust_)
                adjust();

#ifdef USING_GMP
            assert(*this == mpz_class{lhs.template cast<mpz_class>() * rhs.template cast<mpz_class>()});
#endif

            return *this;
        }

        CONSTEXPR Integer& operator+=(Integer const& other)
        {
            if (isNan() || other.isNan())
                setNan();
            else if (isInfinity() || other.isInfinity())
            {

                if ((isPositive() && other.isNegative())
                    || (isNegative() && other.isPositive()))
                    setNan();
                else
                {
                    if (other.isInfinity())
                        isPositive_ = other.isPositive_;

                    setInfinity();
                }
            }
            else
            {
                if ((isPositive() && other.isPositive())
                    || (isNegative() && other.isNegative()))
                {
                    T carry{0};
                    auto const& a{bits_};
                    auto const& b{other.bits_};
                    size_t const n{std::max(a.size(), b.size())};
                    std::vector<T> result;

                    for (size_t i{0}; i < n; ++i)
                    {
                        auto const bit_a{(i < a.size()) ? a[a.size() - 1 - i] : T{0}};
                        auto const bit_b{(i < b.size()) ? b[b.size() - 1 - i] : T{0}};
                        auto const sum{static_cast<T>(bit_a + bit_b + carry)};

                        carry = (sum < bit_a || sum < bit_b);

                        result.emplace_back(sum);
                    }

                    if (carry)
                        result.emplace_back(T{1});

                    std::reverse(result.begin(), result.end());

                    bits_ = result;
                }
                else
                {
                    auto otherBits{other.bits_};

                    if (isPositive())
                    {
                        if (*this < -other)
                        {
                            isPositive_ = false;
                            otherBits = bits_;
                            bits_ = other.bits_;
                        }
                    }
                    else
                    {
                        if (*this > -other)
                        {
                            isPositive_ = true;
                            otherBits = bits_;
                            bits_ = other.bits_;
                        }
                    }

                    T borrow{0};
                    auto const& a{bits_};
                    auto const& b{otherBits};
                    size_t const n{std::max(a.size(), b.size())};
                    std::vector<T> result;

                    for (size_t i{0}; i < n; ++i)
                    {
                        auto const bit_a{(i < a.size()) ? a[a.size() - 1 - i] : T{0}};
                        auto const bit_b{(i < b.size()) ? b[b.size() - 1 - i] : T{0}};
                        auto const bit_result{static_cast<T>(bit_a - bit_b - borrow)};

                        borrow = (bit_result  > bit_a);

                        result.emplace_back(bit_result);
                    }

                    std::reverse(result.begin(), result.end());

                    bits_ = result;
                }
            }

            if (autoAdjust_)
                adjust();

            return *this;
        }

        CONSTEXPR Integer& operator-=(Integer const& other)
        {
            auto const lhs(*this);
            auto const rhs(other);

            *this += -other;

            assert(lhs == *this + rhs);

            return *this;
        }

        CONSTEXPR Integer& operator/=(Integer const& other)
        {
            auto const lhs(*this);
            auto const rhs(other);

            if (other.isNegative())
            {
                *this = -*this;

                return *this /= -other;
            }

            auto const n(*this);

            if (!other || other.isNan())
                setNan();
            else if (other.isInfinity())
                *this = 0;
            else
            {
                if (abs() < other.abs())
                    *this = 0;
                else if (isPositive_ && other.isPositive_)
                {
                    if (this->template fits<longest_type>() && other.template fits<longest_type>())
                        *this = this->template cast<longest_type>() / other.template cast<longest_type>();
                    else if (!(rhs & 1))
                    {
                        auto r(rhs);
                        Integer shift(0);

                        while (!(r & 1))
                        {
                            r >>= 1;
                            ++shift;
                        }

                        *this >>= shift;
                        *this /= r;
                    }
                    else
                        *this = computeQuotientByDivision(*this, other);
                }
                else
                {
                    *this /= -other;
                    *this = -*this;
                }
            }

            assert(abs() <= n.abs());

#ifdef USING_GMP
        assert(*this == mpz_class{lhs.template cast<mpz_class>() / rhs.template cast<mpz_class>()});
#endif

            return *this;
        }

        CONSTEXPR Integer& operator%=(Integer const& other)
        {
            auto const lhs(*this);
            auto const rhs(other);

            if (!other || other.isNan() || other.isInfinity())
                setNan();
            else
            {
                if ((isPositive_ && other.isPositive_) ||
                    (!isPositive_ && !other.isPositive_))
                {
                    if (other == 1)
                        *this = 0;
                    else if (other == 2)
                        *this &= 1;
                    else if (abs().template fits<longest_type>() && other.abs().template fits<longest_type>())
                    {
                        auto const isPositive{isPositive_};

                        *this = abs().template cast<longest_type>() % other.abs().template cast<longest_type>();

                        isPositive_ = isPositive;
                    }
                    else
                    {
                        auto const qr{computeQrByDivision(*this, other)};

                        assert(*this == qr.first * rhs + qr.second);

                        *this = qr.second;
                    }
                }
                else
                {
                    auto const qr{computeQrByDivision(*this, other)};

                    assert(*this == qr.first * rhs + qr.second);

                    *this = qr.second;
                }
            }

#ifdef USING_GMP
            if (lhs > 0 && rhs > 0)
                assert(*this == mpz_class{lhs.template cast<mpz_class>() % rhs.template cast<mpz_class>()});
#endif

            assert(abs() < rhs.abs());

            return *this;
        }

        CONSTEXPR Integer& operator<<=(Integer other)
        {
            assert(other >= 0);

            if (!*this || !other)
                return *this;

            auto const s{static_cast<unsigned short>(sizeof(T) * 8)};
            auto const n(other / s);

            std::vector<T> const v(n.template cast<longest_type>(), T{0});

            bits_.insert(bits_.end(), v.begin(), v.end());

            other -= n * s;

            std::vector<T> bits(bits_.size() + 1, T{0});

            std::copy(bits_.rbegin(), bits_.rend(), bits.rbegin());

            bits_ = bits;

            auto const shift{other.template cast<longest_type>()};

            if (shift)
            {
                for (auto it{bits_.begin() + 1}; it != bits_.end(); ++it)
                {
                    if ((*it >> (sizeof(T) * 8 - shift)))
                        *(it - 1) |= (*it >> (sizeof(T) * 8 - shift));

                    *it <<= shift;
                }
            }

            if (autoAdjust_)
                adjust();

            return *this;
        }

        CONSTEXPR Integer& operator>>=(Integer other)
        {
            assert(other >= 0);

            if (!*this || !other)
                return *this;

            auto const s{static_cast<unsigned short>(sizeof(T) * 8)};
            auto const n(other / s);

            if (bits_.size() < n.template cast<longest_type>())
            {
                bits_ = std::vector<T>{T{0}};

                return *this;
            }

            bits_.resize(bits_.size() - n.template cast<longest_type>());

            other -= n * s;

            auto const shift{other.template cast<longest_type>()};

            if (shift)
            {
                for (auto it{bits_.rbegin()}; it != bits_.rend(); ++it)
                {
                    *it >>= shift;

                    if (it != bits_.rend() - 1 && (*(it + 1) & ((longest_type{1} << shift) - 1)))
                        *it |= (*(it +  1) & ((longest_type{1} << shift) - 1)) << (sizeof(T) * 8 - shift);
                }
            }

            if (bits_.empty())
                bits_.emplace_back(T{0});

            if (autoAdjust_)
                adjust();

            return *this;
        }

        CONSTEXPR bool operator>=(Integer const& other) const
        {
            return operator>(other) || operator==(other);
        }

        CONSTEXPR bool operator>(Integer const& other) const
        {
            if (!isPositive_ && other.isPositive_)
                return false;

            std::vector<T> a(std::max(bits_.size(), other.bits_.size()), T{0});
            std::vector<T> b{a};

            std::copy(bits_.rbegin(), bits_.rend(), a.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), b.rbegin());

            auto const great{a > b};

            return isPositive_ ? great : !great;
        }

        CONSTEXPR bool operator<=(Integer const& other) const
        {
            return operator<(other) || operator==(other);
        }

        CONSTEXPR bool operator<(Integer const& other) const
        {
            if (isPositive_ && !other.isPositive_)
                return false;

            std::vector<T> a(std::max(bits_.size(), other.bits_.size()), T{0});
            std::vector<T> b{a};

            std::copy(bits_.rbegin(), bits_.rend(), a.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), b.rbegin());

            auto const less{a < b};

            return isPositive_ ? less : !less;
        }

        CONSTEXPR bool operator==(Integer const& other) const noexcept
        {
            if (bits_.size() != other.bits_.size())
            {
                if (bits_.size() > other.bits_.size())
                {
                    for (size_t i{0}; i < bits_.size() - other.bits_.size(); ++i)
                    {
                        if (bits_[i])
                            return false;
                    }
                }
                else
                {
                    for (size_t i{0}; i < other.bits_.size() - bits_.size(); ++i)
                    {
                        if (other.bits_[i])
                            return false;
                    }
                }
            }

            bool zero{true};

            auto it1{bits_.rbegin()};
            auto it2{other.bits_.rbegin()};

            for (size_t i{0}; i < std::min(bits_.size(), other.bits_.size()); ++i)
            {
                if (*it1 != *it2)
                    return false;

                if (*it1)
                    zero = false;

                ++it1;
                ++it2;
            }

            if (isPositive_ != other.isPositive_ && !zero)
                return false;

            return true;
        }

        template <typename S>
        CONSTEXPR bool operator==(S const& other) const
        {
            return *this == Integer(other);
        }

        CONSTEXPR bool operator!=(Integer const& other) const
        {
            return !operator==(other);
        }

        template <typename S>
        CONSTEXPR bool operator!=(S const& other) const
        {
            return *this != Integer(other);
        }

        CONSTEXPR Integer operator-() const
        {
            auto x(*this);

            x.isPositive_ = !x.isPositive_;

            return x;
        }

        CONSTEXPR Integer operator~() const
        {
            auto x(*this);

            x.invert();

            return x;
        }

        CONSTEXPR operator bool() const noexcept
        {
            return !!*this;
        }

        CONSTEXPR bool operator!() const noexcept
        {
            for (auto const& b : bits_)
            {
                if (b)
                    return false;
            }

            return true;
        }

        CONSTEXPR Integer& operator--()
        {
            return *this -= 1;
        }

        CONSTEXPR Integer operator--(int)
        {
            auto x(*this);

            operator--();

            return x;
        }

        CONSTEXPR Integer& operator++()
        {
            return *this += 1;
        }

        CONSTEXPR Integer operator++(int)
        {
            auto x(*this);

            operator++();

            return x;
        }

        template <typename S>
        CONSTEXPR Integer& operator+=(S const& other)
        {
            return *this += Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator-=(S const& other)
        {
            return *this -= Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator/=(S const& other)
        {
            return *this /= Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator*=(S const& other)
        {
            return *this *= Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator%=(S const& other)
        {
            return *this %= Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator>>=(S const& other)
        {
            return *this >>= Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator<<=(S const& other)
        {
            return *this <<= Integer(other);
        }

        CONSTEXPR Integer& operator&=(Integer const& other)
        {
            std::vector<T> v1(std::max(bits_.size(), other.bits_.size()), 0);
            std::vector<T> v2{v1};
            std::vector<T> result{v1};

            std::copy(bits_.rbegin(), bits_.rend(), v1.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), v2.rbegin());

            std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T const& a, T const& b) { return a & b; });

            bits_ = result;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        CONSTEXPR Integer& operator|=(Integer const& other)
        {
            std::vector<T> v1(std::max(bits_.size(), other.bits_.size()), 0);
            std::vector<T> v2{v1};
            std::vector<T> result{v1};

            std::copy(bits_.rbegin(), bits_.rend(), v1.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), v2.rbegin());

            std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T const& a, T const& b) { return a | b; });

            bits_ = result;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator|=(S const& other)
        {
            return *this |= Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator&=(S const& other)
        {
            return *this &= Integer(other);
        }

        template <typename S>
        CONSTEXPR Integer& operator^=(S const& other)
        {
            return *this ^= Integer(other);
        }

        CONSTEXPR Integer& operator^=(Integer const& other)
        {
            std::vector<T> v1(std::max(bits_.size(), other.bits_.size()), 0);
            std::vector<T> v2{v1};
            std::vector<T> result{v1};

            std::copy(bits_.rbegin(), bits_.rend(), v1.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), v2.rbegin());

            std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T const& a, T const& b) { return a ^ b; });

            bits_ = result;

            if (autoAdjust_)
                adjust();

            return *this;
        }

        template <typename S>
        CONSTEXPR Integer& operator=(S const& other)
        {
            return *this = Integer(other);
        }

        CONSTEXPR std::string toString(size_t base = 10) const
        {
            assert(2 <= base && base <= 62);

            std::string s;

            if (isNan_)
            {
                if (!isPositive_)
                    s = '-' + s;

                s += "nan";

                return s;
            }
            else if (isInfinity_)
            {
                if (!isPositive_)
                    s = '-' + s;

                s += "inf";

                return s;
            }

            if (base == 2)
            {
    #if __cplusplus >= 202002L
                switch (sizeof(T))
                {
                case 1: //unsigned char
                    for (auto const& b : bits_)
                        s += std::format("{:08b}", b);
                    break;

                case 2: //unsigned short
                    for (auto const& b : bits_)
                        s += std::format("{:016b}", b);
                    break;

                case 4: //unsigned int, unsigned long
                    for (auto const& b : bits_)
                        s += std::format("{:032b}", b);
                    break;

                case 8: //unsigned long long
                    for (auto const& b : bits_)
                        s += std::format("{:064b}", b);
                    break;

                case 16:
                    for (auto const& b : bits_)
                        s += std::format("{:0128b}", b);
                    break;
                }
    #else
                for (auto it{bits_.rbegin()}; it != bits_.rend(); ++it)
                {
                    auto b{*it};

                    for (size_t i{0}; i < sizeof(T) * 8; ++i)
                    {
                        s = (b % 2 ? '1' : '0') + s;
                        b /= 2;
                    }
                }
    #endif

                s = "0b" + s;
            }
            else if (base == 8)
            {
    #if __cplusplus >= 202002L
                if (bits_.size() == 1)
                    s = std::format("{:o}", bits_.back());
                else
    #else
                {
                    auto number(abs());

                    if (!number)
                        s = "0";

                    while (number)
                    {
                        auto const tmp(number % 8);
                        s = std::to_string(tmp.template cast<short>()) + s;
                        number /= 8;
                    }
                }
    #endif

                    s = "0o" + s;
            }
            else if (base == 10)
            {
                auto number(abs());

                if (bits_.size() == 1)
                    s = std::to_string(bits_.back());
                else
                {
                    if (!number)
                        s = "0";

                    auto const n{static_cast<T>(std::log10(static_cast<T>(~T{0})))};
                    T const b(pow(T{10}, n));

                    while (number)
                    {
                        auto const tmp(number % b);
                        std::ostringstream oss;
                        oss << std::setw(n) << std::setfill('0') << tmp.template cast<longest_type>();
                        s = oss.str() + s;
                        number /= b;
                    }

                    size_t i{0};

                    while (s[i] == '0' && i != s.size())
                        ++i;

                    if (i == s.size())
                        i = 0;

                    s = s.substr(i);
                }
            }
            else if (2 < base && base < 16)
            {
                auto number(abs());

                if (bits_.size() == 1)
                    s = std::to_string(bits_.back());
                else
                {
                    if (!number)
                        s = "0";

                    while (number)
                    {
                        auto const tmp(number % static_cast<unsigned char>(base));;
                        s = std::to_string(tmp.template cast<short>()) + s;
                        number /= static_cast<unsigned char>(base);
                    }
                }
            }
            else if (base == 16)
            {
                auto number(abs());
    #if __cplusplus >= 202002L
                switch (sizeof(T))
                {
                case 1: //unsigned char
                    for (auto const& b : bits_)
                        s += std::format("{:02x}", b);
                    break;

                case 2: //unsigned short
                    for (auto const& b : bits_)
                        s += std::format("{:04x}", b);
                    break;

                case 4: //unsigned int, unsigned long
                    for (auto const& b : bits_)
                        s += std::format("{:08x}", b);
                    break;

                case 8: //unsigned long long
                    for (auto const& b : bits_)
                        s += std::format("{:016x}", b);
                    break;

                case 16:
                    for (auto const& b : bits_)
                        s += std::format("{:032x}", b);
                    break;
                }

                if (bits_.size() == 1)
                    s = std::format("{:x}", bits_.back());
                else
    #else
                {
                    if (!number)
                        s = "0";

                    while (number)
                    {
                        auto const tmp(number % 16);
                        if (number < 10)
                            s = std::to_string(tmp.template cast<short>()) + s;
                        else
                            s = (char)('a' + tmp.template cast<short>() - 10) + s;
                        number /= 16;
                    }
                }
    #endif
                    s = "0x" + s;
            }
            else if (base <= 62)
            {
                auto number(abs());

                if (!number)
                    s = "0";

                while (number)
                {
                    auto const tmp(number % 62);
                    if (number < 10)
                        s = std::to_string(tmp.template cast<short>()) + s;
                    else if (number - 10 < 26)
                        s = (char)('a' + tmp.template cast<short>() - 10) + s;
                    else
                        s = (char)('A' + tmp.template cast<short>() - 36) + s;
                    number /= 62;
                }
            }

            if (!isPositive_)
                s = '-' + s;

            return s;
        }

        CONSTEXPR operator char() const noexcept
        {
            return cast<char>();
        }

        CONSTEXPR operator unsigned char() const noexcept
        {
            return cast<unsigned char>();
        }

        CONSTEXPR operator short() const noexcept
        {
            return cast<short>();
        }

        template <typename S>
        CONSTEXPR S cast() const noexcept
        {
            S n{0};

            size_t const iMax{std::min(std::max(longest_type{1}, longest_type{sizeof(S) / sizeof(T)}), longest_type{bits_.size()})};
            auto it{bits_.rbegin() + iMax - 1};

            for (size_t i{0}; i < iMax; ++i)
            {
                n += *it;

                if (i != iMax - 1)
                    n <<= sizeof(T) * 8;

                --it;
            }

            if (!isPositive_)
                n = -n;

            return n;
        }

        CONSTEXPR operator unsigned short() const noexcept
        {
            return cast<unsigned short>();
        }

        CONSTEXPR operator int() const noexcept
        {
            return cast<int>();
        }

        CONSTEXPR operator unsigned int() const noexcept
        {
            return cast<unsigned int>();
        }

        CONSTEXPR operator long() const noexcept
        {
            return cast<long>();
        }

        CONSTEXPR operator unsigned long() const noexcept
        {
            return cast<unsigned long>();
        }

        CONSTEXPR operator long long() const noexcept
        {
            return cast<long long>();
        }

        CONSTEXPR operator unsigned long long() const noexcept
        {
            return cast<unsigned long long>();
        }

        CONSTEXPR bool isNan() const noexcept
        {
            return isNan_;
        }

        CONSTEXPR void setNan() noexcept
        {
            isNan_ = true;
            isInfinity_ = false;
            bits_.clear();
        }

        CONSTEXPR bool isInfinity() const noexcept
        {
            return isInfinity_;
        }

        CONSTEXPR void setInfinity() noexcept
        {
            isNan_ = false;
            isInfinity_ = true;
            bits_.clear();
        }

        CONSTEXPR Integer abs() const
        {
            if (isNegative())
                return -*this;

            return *this;
        }

        CONSTEXPR size_t precision() const noexcept
        {
            return bits_.size();
        }

        CONSTEXPR void setPrecision(size_t precision)
        {
            assert(precision);

            std::vector<T> bits(precision, T{0});

            std::copy(bits_.rbegin(), bits_.rbegin() + std::min(bits_.size(), precision), bits.rbegin());

            bits_ = bits;
        }

        template <typename URNG>
        CONSTEXPR void setRandom(URNG& g) noexcept
        {
            isPositive_ = g() % 2;

            for (auto& b : bits_)
            {
                auto const n{g()};

                if (sizeof(T) <= sizeof(n))
                    b = static_cast<T>(n);
                else
                {
                    auto const iMax{sizeof(T) / sizeof(n)};

                    for (size_t i{0}; i < iMax; ++i)
                    {
                        b += g();

                        if (i != iMax - 1)
                            b <<= sizeof(n) * 8;
                    }
                }
            }
        }

        CONSTEXPR int isPrime(size_t reps = 50) const
        {
            if (*this < 2)
                return 0;
            else if (*this == 2)
                return 2;
            else if(!(*this & 1))
                return 0;

            std::vector<size_t> const primes{3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107};

            auto const sqrtLimit(sqrt(*this));
            auto it{primes.begin()};

            while (it != primes.end() && *it <= sqrtLimit)
            {
                if (!(*this % *it))
                    return 0;

                ++it;
            }

            if (it != primes.end())
                return 2;

            auto s(*this - 1);

            while (!(s & 1))
                s >>= 1;

            auto reduction{[] (Integer const& t, Integer const& R, Integer const& n, Integer const& n_) -> Integer
                {
                    auto const m((t % R) * n_ % R);
                    auto const x((t + m * n) / R);

                    if (x < n)
                        return x;
                    else
                        return x - n;
                }
            };

            auto redmulmod{[&reduction] (Integer const& a, Integer b, Integer const& n,
                                                           Integer const& R, Integer const& n_, Integer const& R2modn) -> Integer
                {
                    auto const reda(reduction(a * R2modn, R, n, n_));
                    auto const redb(reduction(b * R2modn, R, n, n_));
                    auto const redc(reduction(reda * redb, R, n, n_));

                    return reduction(redc, R, n, n_);
                }
            };

            auto mulmod{[] (Integer const& a, Integer b, Integer const& m) -> bool//It returns true if number is prime otherwise false {
                {
                    Integer x(0);
                    auto y{a % m};

                    while (b > 0)
                    {
                        if (b & 1)
                            x = (x + y) % m;

                        y <<= 1;
                        y %= m;
                        b >>= 1;
                    }

                    return x % m;
                }
            };

            auto modulo{[&redmulmod] (Integer const& base, Integer e, Integer const& m,
                                      Integer const& R, Integer const& m_, Integer const& R2modm) -> Integer
                {
                    Integer x(1);
                    auto y{base};

                    while (e > 0)
                    {
                        if (e & 1)
                        {
                            auto const x_(x);
                            x = redmulmod(x, y, m, R, m_, R2modm);
                            while (x < 0)
                                x += m;
                            assert(x == (x_ * y) % m);
                        }

                        auto const y_(y);
                        y = redmulmod(y, y, m, R, m_, R2modm);
                        while (y < 0)
                            y += m;
                        assert(y == (y_ * y_) % m);
                        e >>= 1;
                    }

                    return x % m;
                }
            };

            std::random_device rd;
            auto const number(*this - 1);

            auto const& m(*this);

            Integer R(2);

            while (R <= m)
                R <<= 1;

            if (!(m & 1))
                ++R;

            while (!m.isCoprime(R))
            {
                if (!(m & 1))
                    --R;

                R <<= 1;

                if (!(m & 1))
                    ++R;
            }

            Integer R_, m_;

            auto const d(gcdExtended(R, -m, R_, m_));

            assert(R * R_ - m * m_ == d);

            if (d == -1)
            {
                R_ = -R;
                m_ = -m_;
            }

            auto const R2modm((R * R) % m);

            for (size_t i{0}; i < reps; ++i)
            {
                auto a(*this);
                a.setRandom(rd);
                a.setPositive();
                a = a % number + 1;

                auto temp{s};
                auto mod{modulo(a, temp, *this, R, m_, R2modm)};

                while (temp != number && !mod && mod != number)
                {
                    mod = mulmod(mod, mod, *this);
                    temp <<= 1;
                }

                if (mod != number && !(temp & 1))
                    return 0;
            }

            return 1;
        }

        CONSTEXPR void setPositive()
        {
            isPositive_ = true;
        }

        CONSTEXPR void setNegative()
        {
            isPositive_ = false;
        }

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

        CONSTEXPR bool bit(size_t n) const noexcept
        {
            auto it{bits_.rbegin()};

            while (it != bits_.rend() && n > sizeof(T) * 8)
            {
                n -= sizeof(T) * 8;
                ++it;
            }

            if (it == bits_.rend())
                return false;

            return *it & (T{1} << n);
        }

        CONSTEXPR void setBit(size_t n, bool bit)
        {
            auto it{bits_.rbegin()};

            while (n > sizeof(T) * 8)
            {
                n -= sizeof(T) * 8;
                ++it;

                if (it == bits_.rend())
                {
                    std::vector<T> bits(bits_.size() + 1, T{0});

                    std::copy(bits_.rbegin(), bits_.rend(), bits.rbegin());

                    bits_ = bits;

                    it = bits_.rend() - 1;
                }
            }

            if (bit)
                *it |= T{1} << n;
            else
                *it &= ~(T{1} << n);
        }

        CONSTEXPR T bits(size_t n) const noexcept
        {
            if (n >= bits_.size())
                return T{0};

            return bits_[bits_.size() - 1 - n];
        }

        CONSTEXPR void setBits(size_t n, T const& bits)
        {
            if (bits_.size() < n)
            {
                std::vector<T> bits(n, T{0});

                std::copy(bits_.rbegin(), bits_.rend(), bits.rbegin());

                bits_ = bits;
            }

            bits_[bits_.size() - 1 - n] = bits;

            if (autoAdjust)
                adjust();
        }

        template <size_t N>
        CONSTEXPR void setBits(size_t n, std::bitset<N> const& bits)
        {
            for (size_t i{0}; i < bits.size(); ++i)
                setBit(n + i, bits[i]);

            if (autoAdjust)
                adjust();
        }

        CONSTEXPR size_t count() const noexcept
        {
            size_t count{0};

            for (auto b : bits_)
            {
                while (b)
                {
                    if (b & 1)
                        ++count;

                    b >>= 1;
                }
            }

            return count;
        }

        CONSTEXPR size_t number() const noexcept
        {
            size_t number{0};

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

                number += (std::distance(it, bits_.end()) - 1) * sizeof(T) * 8;
            }

            return number;
        }

        CONSTEXPR bool isEven() const noexcept
        {
            if (bits_.empty())
                return false;

            return !(bits_.back() & 1);
        }

        CONSTEXPR bool isOdd() const noexcept
        {
            if (bits_.empty())
                return false;

            return bits_.back() & 1;
        }

        template <typename S>
        CONSTEXPR bool fits() const
        {
            return (*this == this->template cast<S>());
        }

        CONSTEXPR Integer sign() const
        {
            if (*this < 0)
                return Integer(-1);

            return Integer(1);
        }

        CONSTEXPR void setSign(Integer const& other) noexcept
        {
            isPositive_ = other.isPositive_;
        }

        CONSTEXPR Integer previousPrime() const
        {
            if (isNan())
                return *this;

            if (isInfinity() || *this < 2)
                return nan();

            if (*this == 2)
                return Integer(2);
            else if (*this == 3)
                return Integer(2);

            auto n(*this - 2);

            while (!n.isPrime())
                n -= 2;

            return n;
        }

        CONSTEXPR Integer nextPrime() const
        {
            if (isNan())
                return *this;

            if (*this < 2)
                return Integer(2);
            else if (*this == 2)
                return Integer(3);
            else if (isInfinity())
                return nan();

            auto n(*this + 2);

            while (!n.isPrime())
                n += 2;

            return n;
        }

        CONSTEXPR size_t size() const noexcept
        {
            return bits_.size();
        }

        CONSTEXPR void adjust()
        {
            if (bits_.empty())
                return;

            auto it{bits_.begin()};

            while (!*it && it != bits_.end())
                ++it;

            if (it == bits_.end())
                it = bits_.end() - 1;

            if (it != bits_.begin())
                bits_ = std::vector<T>{it, bits_.end()};
        }

        CONSTEXPR bool isCoprime(Integer const& other) const noexcept
        {
            return gcd(*this, other) == 1;
        }

        CONSTEXPR bool autoAdjust() const noexcept
        {
            return autoAdjust_;
        }

        CONSTEXPR void setAutoAdjust(bool autoAdjust) noexcept
        {
            autoAdjust_ = autoAdjust;
        }

    private:
        bool isPositive_{true};
        std::vector<T> bits_;
        bool isNan_{false};
        bool isInfinity_{false};
        bool autoAdjust_{true};
};

using Integerc = Integer<unsigned char>;
using Integers = Integer<unsigned short>;
using Integeri = Integer<unsigned int>;
using Integerl = Integer<unsigned long>;
using Integerll = Integer<unsigned long long>;
using Integer8 = Integer<uint8_t>;
using Integer16 = Integer<uint16_t>;
using Integer32 = Integer<uint32_t>;
using Integer64 = Integer<uint64_t>;

#ifdef USING_GMP
template <>
template <>
mpz_class Integer<unsigned char>::cast<mpz_class>() const noexcept
{
    auto s{toString(2)};
    s.replace(s.find("0b"), 2, "");

    return mpz_class{s, 2};
}

template <>
template <>
mpz_class Integer<unsigned short>::cast<mpz_class>() const noexcept
{
    auto s{toString(2)};
    s.replace(s.find("0b"), 2, "");

    return mpz_class{s, 2};
}

template <>
template <>
mpz_class Integer<unsigned int>::cast<mpz_class>() const noexcept
{
    auto s{toString(2)};
    s.replace(s.find("0b"), 2, "");

    return mpz_class{s, 2};
}

template <>
template <>
mpz_class Integer<unsigned long>::cast<mpz_class>() const noexcept
{
    auto s{toString(2)};
    s.replace(s.find("0b"), 2, "");

    return mpz_class{s, 2};
}

template <>
template <>
mpz_class Integer<unsigned long long>::cast<mpz_class>() const noexcept
{
    auto s{toString(2)};
    s.replace(s.find("0b"), 2, "");

    return mpz_class{s, 2};
}
#endif

template <>
template <>
std::string Integer<unsigned char>::cast<std::string>() const noexcept
{
    return toString();
}

template <>
template <>
std::string Integer<unsigned short>::cast<std::string>() const noexcept
{
    return toString();
}

template <>
template <>
std::string Integer<unsigned int>::cast<std::string>() const noexcept
{
    return toString();
}

template <>
template <>
std::string Integer<unsigned long>::cast<std::string>() const noexcept
{
    return toString();
}

template <>
template <>
std::string Integer<unsigned long long>::cast<std::string>() const noexcept
{
    return toString();
}

template <typename T>
CONSTEXPR Integer<T> operator*(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs *= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator+(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs += rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator-(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs -= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator/(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs /= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator%(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs %= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator&(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs &= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator|(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs |= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator^(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs ^= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator<<(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs <<= rhs;
}

template <typename T>
CONSTEXPR Integer<T> operator>>(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs >>= rhs;
}

template <typename T, typename S>
CONSTEXPR bool operator>(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator>(Integer<T>(rhs));
}

template <typename T, typename S>
CONSTEXPR bool operator>(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator<(Integer<T>(lhs));
}

template <typename T, typename S>
CONSTEXPR bool operator>=(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator>=(Integer<T>(rhs));
}

template <typename T, typename S>
CONSTEXPR bool operator>=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator<=(Integer<T>(lhs));
}

template <typename T, typename S>
CONSTEXPR bool operator<(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator<(Integer<T>(rhs));
}

template <typename T, typename S>
CONSTEXPR bool operator<(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator>(Integer<T>(lhs));
}

template <typename T, typename S>
CONSTEXPR bool operator<=(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator<=(Integer<T>(rhs));
}

template <typename T, typename S>
CONSTEXPR bool operator<=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator>=(Integer<T>(lhs));
}

template <typename T, typename S>
CONSTEXPR bool operator==(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator==(Integer<T>(lhs));
}

template <typename T, typename S>
CONSTEXPR bool operator!=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator!=(Integer<T>(lhs));
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator+(Integer<T> lhs, S const& rhs)
{
    return lhs += Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator+(S const& lhs, Integer<T> rhs)
{
    return rhs += Integer<T>(lhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator-(Integer<T> lhs, S const& rhs)
{
    return lhs -= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator-(S const& lhs, Integer<T> rhs)
{
    return -(rhs -= Integer<T>(lhs));
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator/(Integer<T> lhs, S const& rhs)
{
    return lhs /= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator/(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) /= rhs;
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator*(Integer<T> lhs, S const& rhs)
{
    return lhs *= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator*(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) *= rhs;
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator%(Integer<T> lhs, S const& rhs)
{
    return lhs %= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator%(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) %= rhs;
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator<<(Integer<T> lhs, S const& rhs)
{
    return lhs <<= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator<<(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) <<= rhs;
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator>>(Integer<T> lhs, S const& rhs)
{
    return lhs >>= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator>>(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) >>= rhs;
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator&(Integer<T> lhs, S const& rhs)
{
    return lhs &= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator&(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) &= rhs;
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator|(Integer<T> lhs, S const& rhs)
{
    return lhs |= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator|(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) |= rhs;
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator^(Integer<T> lhs, S const& rhs)
{
    return lhs ^= Integer<T>(rhs);
}

template <typename T, typename S>
CONSTEXPR Integer<T> operator^(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) ^= rhs;
}

template <typename T>
inline CONSTEXPR std::ostream& operator<<(std::ostream& os, Integer<T> const& n)
{
    return os << n.toString();
}

template <typename T>
CONSTEXPR Integer<T> gcd(Integer<T> const& a, Integer<T> const& b)
{
    if (a.isNan() || b.isNan() || a.isInfinity() || b.isInfinity())
        return Integer<T>::nan();

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
        return gcd(a, b >> 1);
    else if (a.isEven() && b.isOdd())
        return gcd(a >> 1, b);
    else //if (a.isOdd() && b.isOdd())
        return gcd((a - b) >> 1, b);
}

template <typename T, typename S>
CONSTEXPR Integer<T> gcd(Integer<T> const& a, S const& b)
{
    return gcd(a, Integer<T>(b));
}

template <typename T, typename S>
CONSTEXPR Integer<T> gcd(S const& a, Integer<T> const& b)
{
    return gcd(Integer<T>(a), b);
}

template <typename T>
CONSTEXPR Integer<T> lcm(Integer<T> const& a, Integer<T> const& b)
{
    return (a * b).abs() / gcd(a, b);
}

template <typename T, typename S>
CONSTEXPR Integer<T> lcm(Integer<T> const& a, S const& b)
{
    return lcm(a, Integer<T>(b));
}

template <typename T, typename S>
CONSTEXPR Integer<T> lcm(S const& a, Integer<T> const& b)
{
    return lcm(Integer<T>(a), b);
}

template <typename T>
CONSTEXPR Integer<T> gcdExtended(Integer<T> a, Integer<T> b, Integer<T>& u, Integer<T>& v)
{
    if (!a && !b)
        return Integer<T>(0);

    Integer<T> r1(a), u1(1), v1(0);
    Integer<T> r2(b), u2(0), v2(1);
    Integer<T> q, r_temp, u_temp, v_temp;

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

template <typename T, typename S1, typename S2>
CONSTEXPR Integer<T> gcdExtended(S1 const& a, S2 const& b, Integer<T>& u, Integer<T>& v)
{
    return gcdExtended(Integer<T>(a), Integer<T>(b), u, v);
}

template <typename T>
CONSTEXPR Integer<T> factorial(Integer<T> const& n)
{
    if (n.isNan() || n.isInfinity())
        return n;

    assert(n >= 0);

    if (n == 0)
        return Integer<T>(1);

    return n * factorial(n - 1);
}

template <typename T>
CONSTEXPR Integer<T> doubleFactorial(Integer<T> const& n)
{
    return factorial(factorial(n));
}

template <typename T>
CONSTEXPR Integer<T> multiFactorial(Integer<T> const& n, Integer<T> const& m)
{
    return pow(factorial(n), m);
}

template <typename T, typename S>
CONSTEXPR Integer<T> multiFactorial(Integer<T> const& n, S const& m)
{
    return multiFactorial(n, Integer<T>(m));
}

template <typename T, typename S>
CONSTEXPR Integer<T> multiFactorial(S const& n, Integer<T> const& m)
{
    return multiFactorial(Integer<T>(n), m);
}

template <typename T>
CONSTEXPR Integer<T> pow(Integer<T> base, Integer<T> exp)
{
    assert(exp >= 0);

    if (base.isInfinity() || base.isNan())
        return base;

    if (exp.isNan())
        return exp;

    if (exp.isInfinity())
        return exp;

    if (base < 0)
    {
        auto n(pow(base.abs(), exp));

        if (exp & 1)
            n = -n;

        return n;
    }

    if (base == 2)
        return Integer<T>(1) << exp;

    Integer<T> result(1);

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

template <typename T, typename S>
CONSTEXPR Integer<T> pow(Integer<T> const& base, S const& exp)
{
    return pow(base, Integer<T>(exp));
}

template <typename T, typename S>
CONSTEXPR Integer<T> pow(S const& base, Integer<T> const& exp)
{
    return pow(Integer<T>(base), exp);
}

template <typename T>
CONSTEXPR Integer<T> powm(Integer<T> base, Integer<T> exp, Integer<T> const& mod)
{
    assert(exp >= 0);

    Integer<T> result(1);

    auto base_mod(base % mod);

    while (exp > 0)
    {
        if (exp % 2 == 1)
        {
            result *= base_mod;
            result %= mod;
        }

        base_mod *= base_mod;
        base_mod %= mod;

        exp /= 2;
    }

    return result;
}

template <typename T, typename S, typename U>
CONSTEXPR Integer<T> powm(Integer<T> const& base, S const& exp, U const& mod)
{
    return powm(base, Integer<T>(exp), Integer<T>(mod));
}

template <typename T>
CONSTEXPR Integer<T> abs(Integer<T> const& n)
{
    return n.abs();
}

template <typename T>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQr(Integer<T> const& dividend, Integer<T> const& divisor)
{
    if (!divisor)
        return {Integer<T>::nan(), Integer<T>::nan()};
    else if (!dividend)
        return {Integer<T>{0}, Integer<T>{0}};
    else if (dividend.isNan())
        return {dividend, dividend};
    else if (divisor.isNan())
        return {divisor, divisor};
    else if (dividend.isInfinity() || divisor.isInfinity())
        return {Integer<T>::nan(), Integer<T>::nan()};
    else if (divisor.abs() > dividend.abs())
        return {Integer<T>{0}, dividend};
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

    Integer<T> start(1);
    auto end(dividend);

    while (start <= end)
    {
        auto mid(end + start);
        mid >>= 1;

        auto n(dividend - divisor * mid);

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

            return {mid, n};
        }
    }

    return {Integer<T>(0), dividend};
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQr(Integer<T> const& dividend, S const& divisor)
{
    return computeQr(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQr(S const& dividend, Integer<T> const& divisor)
{
    return computeQr(Integer<T>(dividend), divisor);
}

template <typename T>
CONSTEXPR Integer<T> computeQuotient(Integer<T> const& dividend, Integer<T> const& divisor)
{
    if (!divisor)
        return Integer<T>::nan();
    else if (!dividend)
        return Integer<T>{0};
    else if (dividend.isNan())
        return dividend;
    else if (divisor.isNan())
        return divisor;
    else if (dividend.isInfinity() || divisor.isInfinity())
        return Integer<T>::nan();
    else if (divisor.abs() > dividend.abs())
        return Integer<T>{0};
    else if (dividend < 0 && divisor < 0)
        return computeQuotient(-dividend, -divisor);
    else if (dividend > 0 && divisor < 0)
        return -computeQuotient(dividend, -divisor);
    else if (dividend < 0 && divisor > 0)
        return -computeQuotient(-dividend, divisor);

    Integer<T> start(1);
    auto end(dividend);

    while (start <= end)
    {
        auto mid(end + start);
        mid >>= 1;

        auto n(dividend - divisor * mid);

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

    return Integer<T>(0);
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQuotient(Integer<T> const& dividend, S const& divisor)
{
    return computeQuotient(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQuotient(S const& dividend, Integer<T> const& divisor)
{
    return computeQuotient(Integer<T>(dividend), divisor);
}

template <typename T, typename S>
CONSTEXPR Integer<T> fibonacci(Integer<T> n)
{
    assert(n >= 0);

    if (!n)
        return 0;
    else if (n == 1)
        return 1;

    n -= 1;

    Integer<T> fn_2(0);
    Integer<T> fn_1(1);
    Integer<T> fn;

    while (n)
    {
        fn = fn_1 + fn_2;
        fn_2 = fn_1;
        fn_1 = fn;
        --n;
    }

    return fn;
}

template <typename T>
CONSTEXPR Integer<T> primorial(Integer<T> n)
{
    Integer<T> result(1);
    Integer<T> number(2);

    while (number <= n)
    {
        result *= number;
        number = number.nextPrime();
    }

    return result;
}

template <typename T>
CONSTEXPR int jacobi(Integer<T> const& a, Integer<T> n)
{
    assert(n > 0 && n.isOdd());

    int result(1);
    Integer<T> prime(2);

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

template <typename T, typename S>
CONSTEXPR int jacobi(Integer<T> const& a, S const& n)
{
    return jacobi(a, Integer<T>(n));
}

template <typename T, typename S>
CONSTEXPR int jacobi(S const& a, Integer<T> const& n)
{
    return jacobi(Integer<T>(a), n);
}

template <typename T>
CONSTEXPR int legendre(Integer<T> const& a, Integer<T> const& p)
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

template <typename T, typename S>
CONSTEXPR Integer<T> legendre(Integer<T> const& a, S const& p)
{
    return legendre(a, Integer<T>(p));
}

template <typename T, typename S>
CONSTEXPR Integer<T> legendre(S const& a, Integer<T> const& p)
{
    return legendre(Integer<T>(a), p);
}

template <typename T>
CONSTEXPR int kronecker(Integer<T> const& a, Integer<T> const& b)
{
    if (a == b)
        return 1;
    else
        return 0;
}

template <typename T, typename S>
CONSTEXPR int kronecker(Integer<T> const& a, S const& b)
{
    return kronecker(a, Integer<T>(b));
}

template <typename T, typename S>
CONSTEXPR int kronecker(S const& a, Integer<T> const& b)
{
    return kronecker(Integer<T>(a), b);
}

template <typename T>
CONSTEXPR Integer<T> binomial(Integer<T> const& n, Integer<T> const& k)
{
    assert(n >= 0 && k >= 0);

    return factorial(n) / (factorial(k) * factorial(n - k));
}

template <typename T, typename S>
CONSTEXPR Integer<T> binomial(Integer<T> const& n, S const& k)
{
    return binomial(n, Integer<T>(k));
}

template <typename T, typename S>
CONSTEXPR Integer<T> binomial(S const& n, Integer<T> const& k)
{
    return binomial(Integer<T>(n), k);
}

template <typename T>
CONSTEXPR Integer<T> sqrt(Integer<T> const& n)
{
    if (n < 0)
        return Integer<T>::nan();
    else if (!n || n == 1)
        return n;

    Integer<T> lo(1), hi(n);
    Integer<T> res(1);

    while (lo <= hi)
    {
        auto mid(lo + hi);
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

template <typename T>
CONSTEXPR Integer<T> computeQuotientBinary(Integer<T> dividend, Integer<T> const& divisor)
{
    if (!divisor)
        return Integer<T>::nan();
    else if (!dividend)
        return Integer<T>{0};
    else if (dividend.isNan())
        return dividend;
    else if (divisor.isNan())
        return divisor;
    else if (dividend.isInfinity() || divisor.isInfinity())
        return Integer<T>::nan();
    else if (divisor.abs() > dividend.abs())
        return Integer<T>{0};
    else if (dividend < 0 && divisor > 0)
        return -computeQuotientBinary(-dividend, divisor);
    else if (dividend > 0 && divisor < 0)
        return -computeQuotientBinary(dividend, -divisor);
    else if (dividend < 0 && divisor < 0)
        return computeQuotientBinary(-dividend, -divisor);

    Integer<T> quotient(0);
    auto tempDivisor(divisor);
    Integer<T> bit(1);

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

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQuotientBinary(Integer<T> const& dividend, S const& divisor)
{
    return computeQuotientBinary(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQuotientBinary(S const& dividend, Integer<T> const& divisor)
{
    return computeQuotientBinary(Integer<T>(dividend), divisor);
}

template <typename T>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQrBinary(Integer<T> const& dividend, Integer<T> const& divisor)
{
    std::pair<Integer<T>, Integer<T> > qr{computeQuotientBinary(dividend, divisor), Integer<T>(0)};
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

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQrBinary(Integer<T> const& dividend, S const& divisor)
{
    return computeQrBinary(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQrBinary(S const& dividend, Integer<T> const& divisor)
{
    return computeQrBinary(Integer<T>(dividend), divisor);
}

template <typename T>
CONSTEXPR Integer<T> computeQuotientByDivision(Integer<T> dividend, Integer<T> const& divisor)
{
    if (!divisor)
        return Integer<T>::nan();
    else if (!dividend)
        return Integer<T>{0};
    else if (dividend.isNan())
        return dividend;
    else if (divisor.isNan())
        return divisor;
    else if (dividend.isInfinity() || divisor.isInfinity())
        return Integer<T>::nan();
    else if (divisor.abs() > dividend.abs())
        return Integer<T>{0};

    std::pair<Integer<T>, Integer<T> > qr;
    qr = computeQrByDivision(dividend, divisor);

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

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQuotientByDivision(Integer<T> const& dividend, S const& divisor)
{
    return computeQuotientByDivision(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQuotientByDivision(S const& dividend, Integer<T> const& divisor)
{
    return computeQuotientByDivision(Integer<T>(dividend), divisor);
}

template <typename T>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQrByDivision(Integer<T> const& dividend, Integer<T> const& divisor)
{
    if (!divisor)
        return {Integer<T>::nan(), Integer<T>::nan()};
    else if (!dividend)
        return {Integer<T>{0}, Integer<T>{0}};
    else if (dividend.isNan())
        return {dividend, dividend};
    else if (divisor.isNan())
        return {divisor, divisor};
    else if (dividend.isInfinity() || divisor.isInfinity())
        return {Integer<T>::nan(), Integer<T>::nan()};
    else if (divisor.abs() > dividend.abs())
        return {Integer<T>{0}, dividend};
    else if (dividend < 0 && divisor < 0)
    {
        auto qr{computeQrByDivision(-dividend, -divisor)};

        if (qr.second)
        {
            ++qr.first;
            qr.second += divisor;
        }

        return qr;
    }
    else if (dividend > 0 && divisor < 0)
    {
        auto qr{computeQrByDivision(dividend, -divisor)};

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
        auto qr{computeQrByDivision(-dividend, divisor)};

        qr.first = -qr.first;

        if (qr.second)
        {
            --qr.first;
            qr.second = divisor - qr.second;
        }

        return qr;
    }

    Integer<T> mask;
    mask.setPrecision(sizeof(longest_type) / sizeof(T));
    mask = ~mask;
    mask >>= sizeof(longest_type) * 4;

    std::vector<longest_type> divisors;
    auto divisorTmp(divisor);

    for (size_t i{0}; i < divisor.number() / (sizeof(longest_type) * 4) + 1; ++i)
    {
        divisors.emplace_back((divisorTmp & mask).template cast<longest_type>());
        divisorTmp >>= sizeof(longest_type) * 4;
    }

    while (!divisors.back())
        divisors.pop_back();

    std::reverse(divisors.begin(), divisors.end());

    Integer<T> mask1;
    mask1.setPrecision(dividend.size());
    mask1 = ~mask1;

    mask1 >>= (sizeof(longest_type) * 8) * (divisor.number() / (sizeof(longest_type) * 8))
              + (divisor.number() % (sizeof(longest_type) * 8) ? sizeof(longest_type) * 8 : 0);

    mask1.setPrecision(dividend.size());

    auto mask2(~mask1);
    Integer<T> quotient, remainder;
    auto d(dividend);

    Integer<T> maskLongestType2;
    maskLongestType2.setPrecision(sizeof(longest_type) / sizeof(T));
    maskLongestType2 = ~maskLongestType2;
    maskLongestType2 >>= sizeof(longest_type) * 4;

    while (d >= divisor)
    {
        auto dividendTmp(d);

        dividendTmp = d & mask2;

        Integer<T> shift(0);

        for (auto m(mask2); !(m & 1); m >>= 1)
        {
            dividendTmp >>= 1;
            ++shift;
        }

        std::vector<longest_type> dividends;
        auto dTmp(d);

        for (size_t i{0}; i < d.number() / (sizeof(longest_type) * 4) + 1; ++i)
        {
            dividends.emplace_back((dTmp & maskLongestType2).template cast<longest_type>());
            dTmp >>= sizeof(longest_type) * 4;
        }

        while (!dividends.back())
            dividends.pop_back();

        if (dividends.empty())
            dividends.emplace_back(T{0});

        std::vector<longest_type> dividendsTmp;
        dTmp = dividendTmp;

        for (size_t i{0}; i < d.number() / (sizeof(longest_type) * 4) + 1; ++i)
        {
            dividendsTmp.emplace_back((dTmp & maskLongestType2).template cast<longest_type>());
            dTmp >>= sizeof(longest_type) * 4;
        }

        while (!dividendsTmp.back())
            dividendsTmp.pop_back();

        if (dividendsTmp.empty())
            dividendsTmp.emplace_back(T{0});

        std::reverse(dividends.begin(), dividends.end());

        auto q{dividends.front() / divisors.front()};
        size_t j{0};

        if (!q && dividends.size() > 0)
        {
            q = ((dividends[0] << sizeof(longest_type) * 4) | dividends[1]) / divisors.front();

            if (!q && dividends.size() > 1)
                q = ((((dividends[0] << sizeof(longest_type) * 4) | dividends[1]) << sizeof(longest_type) * 4) | dividends[2]) / divisors.front();
        }

        while (dividendTmp < divisor && j + 1 < dividends.size())
        {
            ++j;
            dividendTmp = (dividendTmp << sizeof(longest_type) * 4) | dividends[j];
        }

        while (q * divisor > dividendTmp)
            --q;

        while ((q + 1) * divisor <= dividendTmp)
            ++q;

        quotient <<= sizeof(longest_type) * 4;
        quotient |= q;

        remainder = dividendTmp - q * divisor;

        dTmp = remainder;

        if (dividends.size() >= divisors.size() + j)
        {
            dTmp <<= (dividends.size() - divisors.size() - j) * sizeof(T) * 8;

            Integer<T> m;

            for (size_t i{0}; i < (dividends.size() - divisors.size() - j) * sizeof(T) * 8; ++i)
            {
                m <<= 1;
                m |= 1;
            }

            dTmp |= dividend & m;
        }

        d = dTmp;
    }

    return {quotient, remainder};
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQrByDivision(Integer<T> const& dividend, S const& divisor)
{
    return computeQrByDivision(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
CONSTEXPR std::pair<Integer<T>, Integer<T> > computeQrByDivision(S const& dividend, Integer<T> const& divisor)
{
    return computeQrByDivision(Integer<T>(dividend), divisor);
}

Integerc operator""_zc(char const* str)
{
    return Integerc(str);
}

Integers operator""_zs(char const* str)
{
    return Integers(str);
}

Integeri operator""_zi(char const* str)
{
    return Integeri(str);
}

Integerl operator""_zl(char const* str)
{
    return Integerl(str);
}

Integerll operator""_zll(char const* str)
{
    return Integerll(str);
}

Integer8 operator""_z8(char const* str)
{
    return Integer8(str);
}

Integer16 operator""_z16(char const* str)
{
    return Integer16(str);
}

Integer32 operator""_z32(char const* str)
{
    return Integer32(str);
}

Integer64 operator""_z64(char const* str)
{
    return Integer64(str);
}

Integer<uintmax_t> operator""_z(char const* str)
{
    return Integer<uintmax_t>(str);
}

#endif // INTEGER_H
