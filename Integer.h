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

#if __cplusplus >= 202002L
#include <format>
#endif

#ifdef USING_GMP
#include <gmpxx.h>
#endif

using longest_type = unsigned long long;

template <typename T, typename Enable = void>
class Integer;

template <typename T>
class Integer<T, typename std::enable_if<std::is_unsigned<T>::value>::type>
{
    public:
        constexpr Integer() : bits_{T{0}}
        {
        }

        template <typename S>
        constexpr explicit Integer(S n) : isPositive_{n >= 0}
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

        constexpr explicit Integer(std::vector<T> const& bits, bool isPositive = true) : isPositive_{isPositive}, bits_{bits}
        {
            if (bits_.empty())
                bits_.emplace_back(T{0});
        }

        template <size_t N>
        constexpr explicit Integer(std::bitset<N> const& bits, bool isPositive = true) : isPositive_{isPositive}
        {
            setBits(0, bits);
        }

        constexpr explicit Integer(std::initializer_list<T> const& bits, bool isPositive = true) : isPositive_{isPositive}, bits_{bits}
        {
            if (bits_.empty())
                bits_.emplace_back(T{0});
        }

        template <class InputIt>
        constexpr explicit Integer(InputIt begin, InputIt end, bool isPositive = true) : isPositive_{isPositive}, bits_{begin, end}
        {
            if (bits_.empty())
                bits_.emplace_back(T{0});
        }

#ifdef USING_GMP
        constexpr explicit Integer(mpz_class const& n) : Integer(n.get_str())
        {
        }
#endif

        constexpr explicit Integer(char const* n, size_t base = 10) : Integer(std::string{n}, base)
        {
        }

        constexpr explicit Integer(std::string n, size_t base = 10)
        {
            assert(2 <= base && base <= 62);

            n.erase(std::remove_if(n.begin(), n.end(), isspace), n.end());
            n.erase(std::remove(n.begin(), n.end(), '\''), n.end());

            auto it{n.begin()};

            if (*it == '-')
            {
                isPositive_ = false;
                ++it;
            }

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
                    else if (n.substr(0, 2) == "0b" || n.substr(0, 2) == "0B")
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
                    else if (n.substr(0, 2) == "0o" || n.substr(0, 2) == "0O")
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
                    else if (n.substr(0, 2) == "0x" || n.substr(0, 2) == "0X")
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
        }

        constexpr bool isPositive() const noexcept
        {
            return isPositive_;
        }

        constexpr bool isNegative() const noexcept
        {
            return !isPositive_;
        }

        constexpr auto const& bits() const noexcept
        {
            return bits_;
        }

        constexpr void invert() noexcept
        {
            for (size_t i{0}; i < bits_.size(); ++i)
                bits_[i] = ~bits_[i];
        }

        constexpr Integer& operator*=(Integer const& other)
        {
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
                        auto const a(this->template cast<longest_type>());
                        auto const b(other.template cast<longest_type>());
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

                            auto const z0(x0 * y0);
                            auto z1(x1 * y0 + x0 * y1);
                            z1 <<= m;
                            auto z2(x1 * y1);
                            z2 <<= 2 * m;

                            //xy = z2 * 2^(2 * m) + z1 * 2^m + z0
                            *this = z0 + z1 + z2;
                        }
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
                                  bits_.rbegin() + std::min(n1 / 2, m),
                                  x0.bits_.rbegin());
                        x1.bits_ = std::vector<T>(m, T{0});
                        std::copy(bits_.rbegin() + std::min(n1 / 2, m),
                                  bits_.rbegin() + std::min(bits_.size(), n1),
                                  x1.bits_.rbegin());
                        y0.bits_ = std::vector<T>(m, T{0});
                        std::copy(other.bits_.rbegin(),
                                  other.bits_.rbegin() + std::min(n2 / 2, m),
                                  y0.bits_.rbegin());
                        y1.bits_ = std::vector<T>(m, T{0});
                        std::copy(other.bits_.rbegin() + std::min(n2 / 2, m),
                                  other.bits_.rbegin() + std::min(other.bits_.size(), n2),
                                  y1.bits_.rbegin());

                        auto const z0(x0 * y0);
                        auto const z1(x1 * y0 + x0 * y1);
                        auto const z2(x1 * y1);

                        //xy = z2 * 2^(2 * m) + z1 * 2^m + z0
                        *this = z0;
                        Integer w1, w2;
                        w1.bits_ = std::vector<T>(z1.bits_.size() + m, T{0});
                        std::copy(z1.bits_.rbegin(), z1.bits_.rend(), w1.bits_.rbegin() + m);
                        *this += w1;
                        w2.bits_ = std::vector<T>(z2.bits_.size() + 2 * m, T{0});
                        std::copy(z2.bits_.rbegin(), z2.bits_.rend(), w2.bits_.rbegin() + 2 * m);
                        *this += w2;

                        adjust();
                    }
                }
                else
                {
                    *this *= -other;
                    *this = -*this;
                }
            }

            return *this;
        }

        constexpr Integer& operator+=(Integer const& other)
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

            return *this;
        }

        constexpr Integer& operator-=(Integer const& other)
        {
            auto const n(*this);

            *this += -other;

            assert(n == *this + other);

            return *this;
        }

        constexpr Integer& operator/=(Integer other)
        {
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
                    else
                        *this = computeQuotient(*this, other);
                }
                else
                {
                    *this /= -other;
                    *this = -*this;
                }
            }

            assert(abs() <= n.abs());

            return *this;
        }

        constexpr Integer& operator%=(Integer const& other)
        {
            std::cout << "*this ";
            for (auto const& b : bits_)
                std::cout << b << " ";
            std::cout << std::endl;
            std::cout << "other ";
            for (auto const& b : other.bits_)
                std::cout << b << " ";
            std::cout << std::endl;
            if (!other || other.isNan() || other.isInfinity())
                setNan();
            else
            {
                if ((isPositive_ && other.isPositive_) ||
                    (!isPositive_ && !other.isPositive_))
                {
                    if (abs().template fits<longest_type>() && other.abs().template fits<longest_type>())
                    {
                        auto const isPositive{isPositive_};

                        *this = abs().template cast<longest_type>() % other.abs().template cast<longest_type>();

                        isPositive_ = isPositive;
                    }
                    else
                        *this = computeQr(*this, other).second;
                }
                else
                    *this = computeQr(*this, other).second;
            }

            assert(abs() < other.abs());

            return *this;
        }

        constexpr Integer& operator<<=(Integer other)
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

            adjust();

            return *this;
        }

        constexpr Integer& operator>>=(Integer other)
        {
            assert(other >= 0);

            if (!*this || !other)
                return *this;

            auto const s{static_cast<unsigned short>(sizeof(T) * 8)};
            auto const n(other / s);

            bits_.resize(bits_.size() - n.template cast<longest_type>());

            other -= n * s;

            auto const shift{other.template cast<longest_type>()};

            if (shift)
            {
                for (auto it{bits_.rbegin()}; it != bits_.rend(); ++it)
                {
                    *it >>= shift;

                    if (it != bits_.rend() - 1 && (*(it + 1) & ((1 << shift) - 1)))
                        *it |= (*(it +  1) & ((1 << shift) - 1)) << (sizeof(T) * 8 - shift);
                }
            }

            if (bits_.empty())
                bits_.emplace_back(T{0});

            adjust();

            return *this;
        }

        constexpr bool operator>=(Integer const& other) const
        {
            return operator>(other) || operator==(other);
        }

        constexpr bool operator>(Integer const& other) const
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

        constexpr bool operator<=(Integer const& other) const
        {
            return operator<(other) || operator==(other);
        }

        constexpr bool operator<(Integer const& other) const
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

        constexpr bool operator==(Integer const& other) const noexcept
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
        constexpr bool operator==(S const& other) const
        {
            return *this == Integer(other);
        }

        constexpr bool operator!=(Integer const& other) const
        {
            return !operator==(other);
        }

        template <typename S>
        constexpr bool operator!=(S const& other) const
        {
            return *this != Integer(other);
        }

        constexpr Integer operator-() const
        {
            auto x(*this);

            x.isPositive_ = !x.isPositive_;

            return x;
        }

        constexpr Integer operator~() const
        {
            auto x(*this);

            x.invert();

            return x;
        }

        constexpr operator bool() const noexcept
        {
            return !!*this;
        }

        constexpr bool operator!() const noexcept
        {
            for (auto const& b : bits_)
            {
                if (b)
                    return false;
            }

            return true;
        }

        constexpr Integer& operator--()
        {
            return *this -= 1;
        }

        constexpr Integer operator--(int)
        {
            auto x(*this);

            operator--();

            return x;
        }

        constexpr Integer& operator++()
        {
            return *this += 1;
        }

        constexpr Integer operator++(int)
        {
            auto x(*this);

            operator++();

            return x;
        }

        template <typename S>
        constexpr Integer& operator+=(S const& other)
        {
            return *this += Integer(other);
        }

        template <typename S>
        constexpr Integer& operator-=(S const& other)
        {
            return *this -= Integer(other);
        }

        template <typename S>
        constexpr Integer& operator/=(S const& other)
        {
            return *this /= Integer(other);
        }

        template <typename S>
        constexpr Integer& operator*=(S const& other)
        {
            return *this *= Integer(other);
        }

        template <typename S>
        constexpr Integer& operator%=(S const& other)
        {
            return *this %= Integer(other);
        }

        template <typename S>
        constexpr Integer& operator>>=(S const& other)
        {
            return *this >>= Integer(other);
        }

        template <typename S>
        constexpr Integer& operator<<=(S const& other)
        {
            return *this <<= Integer(other);
        }

        constexpr Integer& operator&=(Integer const& other)
        {
            std::vector<T> v1(std::max(bits_.size(), other.bits_.size()), 0);
            std::vector<T> v2{v1};
            std::vector<T> result{v1};

            std::copy(bits_.rbegin(), bits_.rend(), v1.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), v2.rbegin());

            std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T const& a, T const& b) { return a & b; });

            bits_ = result;

            return *this;
        }

        constexpr Integer& operator|=(Integer const& other)
        {
            std::vector<T> v1(std::max(bits_.size(), other.bits_.size()), 0);
            std::vector<T> v2{v1};
            std::vector<T> result{v1};

            std::copy(bits_.rbegin(), bits_.rend(), v1.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), v2.rbegin());

            std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T const& a, T const& b) { return a | b; });

            bits_ = result;

            return *this;
        }

        template <typename S>
        constexpr Integer& operator|=(S const& other)
        {
            return *this |= Integer(other);
        }

        template <typename S>
        constexpr Integer& operator&=(S const& other)
        {
            return *this &= Integer(other);
        }

        template <typename S>
        constexpr Integer& operator^=(S const& other)
        {
            return *this ^= Integer(other);
        }

        constexpr Integer& operator^=(Integer const& other)
        {
            std::vector<T> v1(std::max(bits_.size(), other.bits_.size()), 0);
            std::vector<T> v2{v1};
            std::vector<T> result{v1};

            std::copy(bits_.rbegin(), bits_.rend(), v1.rbegin());
            std::copy(other.bits_.rbegin(), other.bits_.rend(), v2.rbegin());

            std::transform(v1.begin(), v1.end(), v2.begin(), result.begin(), [](T const& a, T const& b) { return a ^ b; });

            bits_ = result;

            return *this;
        }

        template <typename S>
        constexpr Integer& operator=(S const& other)
        {
            return *this = Integer(other);
        }

        constexpr std::string toString(size_t base = 10) const
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

        constexpr operator char() const noexcept
        {
            return cast<char>();
        }

        constexpr operator unsigned char() const noexcept
        {
            return cast<unsigned char>();
        }

        constexpr operator short() const noexcept
        {
            return cast<short>();
        }

        template <typename S>
        constexpr S cast() const noexcept
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

        constexpr operator unsigned short() const noexcept
        {
            return cast<unsigned short>();
        }

        constexpr operator int() const noexcept
        {
            return cast<int>();
        }

        constexpr operator unsigned int() const noexcept
        {
            return cast<unsigned int>();
        }

        constexpr operator long() const noexcept
        {
            return cast<long>();
        }

        constexpr operator unsigned long() const noexcept
        {
            return cast<unsigned long>();
        }

        constexpr operator long long() const noexcept
        {
            return cast<long long>();
        }

        constexpr operator unsigned long long() const noexcept
        {
            return cast<unsigned long long>();
        }

        constexpr bool isNan() const noexcept
        {
            return isNan_;
        }

        constexpr void setNan() noexcept
        {
            isNan_ = true;
            isInfinity_ = false;
            bits_.clear();
        }

        constexpr bool isInfinity() const noexcept
        {
            return isInfinity_;
        }

        constexpr void setInfinity() noexcept
        {
            isNan_ = false;
            isInfinity_ = true;
            bits_.clear();
        }

        constexpr Integer abs() const
        {
            if (isNegative())
                return -*this;

            return *this;
        }

        constexpr size_t precision() const noexcept
        {
            return bits_.size();
        }

        constexpr void setPrecision(size_t precision)
        {
            assert(precision);

            std::vector<T> bits(precision, T{0});

            std::copy(bits_.rbegin(), bits_.rbegin() + std::min(bits_.size(), precision), bits.rbegin());

            bits_ = bits;
        }

        template <typename URNG>
        constexpr void setRandom(URNG& g) noexcept
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

        constexpr bool isPrime(size_t reps = 50) const
        {
            if (*this < 2)
                return false;

            if (*this != 2 && *this % 2 == 0)
                return false;

            auto s(*this - 1);

            while (s % 2 == 0)
                s /= 2;

            auto mulmod{[] (Integer const& a, Integer b, Integer const& m) -> Integer//It returns true if number is prime otherwise false {
                {
                    Integer x{0};
                    auto y{a % m};

                    while (b > 0)
                    {
                        if (b % 2 == 1)
                            x = (x + y) % m;

                        y = (y * 2) % m;
                        b /= 2;
                    }

                    return x % m;
                }
            };

            auto modulo{[] (Integer const& base, Integer e, Integer const& m) -> Integer
                {
                    Integer x{1};
                    auto y{base};

                    while (e > 0)
                    {
                        if (e % 2 == 1)
                            x = (x * y) % m;

                        y = (y * y) % m;
                        e = e / 2;
                    }

                    return x % m;
                }
            };

            std::random_device rd;
            auto const number(*this - 1);

            for (size_t i{0}; i < reps; ++i)
            {
                auto a(*this);
                a.setRandom(rd);
                a.setPositive();
                a = a % number + 1;

                auto temp{s};
                auto mod{modulo(a, temp, *this)};

                while (temp != number && !mod && mod != number)
                {
                    mod = mulmod(mod, mod, *this);
                    temp *= 2;
                }

                if (mod != number && temp % 2 == 0)
                    return false;
            }

            return true;
        }

        constexpr void setPositive()
        {
            isPositive_ = true;
        }

        constexpr void setNegative()
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

        constexpr bool bit(size_t n) const noexcept
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

        constexpr void setBit(size_t n, bool bit)
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

        constexpr T bits(size_t n) const noexcept
        {
            if (n >= bits_.size())
                return T{0};

            return bits_[bits_.size() - 1 - n];
        }

        constexpr void setBits(size_t n, T const& bits)
        {
            if (bits_.size() < n)
            {
                std::vector<T> bits(n, T{0});

                std::copy(bits_.rbegin(), bits_.rend(), bits.rbegin());

                bits_ = bits;
            }

            bits_[bits_.size() - 1 - n] = bits;
        }

        template <size_t N>
        constexpr void setBits(size_t n, std::bitset<N> const& bits)
        {
            for (size_t i{0}; i < bits.size(); ++i)
                setBit(n + i, bits[i]);
        }

        constexpr size_t count() const noexcept
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

        constexpr size_t number() const noexcept
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

        constexpr bool isEven() const noexcept
        {
            if (bits_.empty())
                return false;

            return !(bits_.back() & 1);
        }

        constexpr bool isOdd() const noexcept
        {
            if (bits_.empty())
                return false;

            return bits_.back() & 1;
        }

        template <typename S>
        constexpr bool fits() const
        {
            return (*this == this->template cast<S>());
        }

        constexpr Integer sign() const
        {
            if (*this < 0)
                return Integer(-1);

            return Integer(1);
        }

        constexpr void setSign(Integer const& other) noexcept
        {
            isPositive_ = other.isPositive_;
        }

        constexpr Integer previousPrime() const
        {
            if (isNan())
                return *this;

            if (isInfinity() || *this < 2)
                return nan();

            if (*this == 2)
                return 2;
            else if (*this == 3)
                return 2;

            auto n(*this - 2);

            while (!n.isPrime())
                n -= 2;
        }

        constexpr Integer nextPrime() const
        {
            if (isNan())
                return *this;

            if (*this < 2)
                return 2;
            else if (*this == 2)
                return 3;
            else if (isInfinity())
                return nan();

            auto n(*this + 2);

            while (!n.isPrime())
                n += 2;
        }

        constexpr size_t size() const noexcept
        {
            return bits_.size();
        }

        constexpr void adjust()
        {
            if (bits_.empty())
                return;

            auto it{bits_.begin()};

            while (!*it && it != bits_.end())
                ++it;

            if (it == bits_.end())
                it = bits_.end() - 1;

            bits_ = std::vector<T>{it, bits_.end()};
        }

#ifdef USING_GMP
        constexpr mpz_class toMpz_class() const
        {
            return mpz_class{toString(2).substr(2), 2};
        }
#endif
    private:
        bool isPositive_{true};
        std::vector<T> bits_;
        bool isNan_{false};
        bool isInfinity_{false};
};

template <typename T>
constexpr Integer<T> operator*(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs *= rhs;
}

template <typename T>
constexpr Integer<T> operator+(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs += rhs;
}

template <typename T>
constexpr Integer<T> operator-(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs -= rhs;
}

template <typename T>
constexpr Integer<T> operator/(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs /= rhs;
}

template <typename T>
constexpr Integer<T> operator%(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs %= rhs;
}

template <typename T>
constexpr Integer<T> operator&(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs &= rhs;
}

template <typename T>
constexpr Integer<T> operator|(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs |= rhs;
}

template <typename T>
constexpr Integer<T> operator^(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs ^= rhs;
}

template <typename T>
constexpr Integer<T> operator<<(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs <<= rhs;
}

template <typename T>
constexpr Integer<T> operator>>(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs >>= rhs;
}

template <typename T, typename S>
constexpr bool operator>(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator>(Integer<T>(rhs));
}

template <typename T, typename S>
constexpr bool operator>(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator<(Integer<T>(lhs));
}

template <typename T, typename S>
constexpr bool operator>=(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator>=(Integer<T>(rhs));
}

template <typename T, typename S>
constexpr bool operator>=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator<=(Integer<T>(lhs));
}

template <typename T, typename S>
constexpr bool operator<(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator<(Integer<T>(rhs));
}

template <typename T, typename S>
constexpr bool operator<(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator>(Integer<T>(lhs));
}

template <typename T, typename S>
constexpr bool operator<=(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator<=(Integer<T>(rhs));
}

template <typename T, typename S>
constexpr bool operator<=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator>=(Integer<T>(lhs));
}

template <typename T, typename S>
constexpr bool operator==(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator==(Integer<T>(lhs));
}

template <typename T, typename S>
constexpr bool operator!=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator!=(Integer<T>(lhs));
}

template <typename T, typename S>
constexpr Integer<T> operator+(Integer<T> lhs, S const& rhs)
{
    return lhs += Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator+(S const& lhs, Integer<T> rhs)
{
    return rhs += Integer<T>(lhs);
}

template <typename T, typename S>
constexpr Integer<T> operator-(Integer<T> lhs, S const& rhs)
{
    return lhs -= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator-(S const& lhs, Integer<T> rhs)
{
    return -(rhs -= Integer<T>(lhs));
}

template <typename T, typename S>
constexpr Integer<T> operator/(Integer<T> lhs, S const& rhs)
{
    return lhs /= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator/(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) /= rhs;
}

template <typename T, typename S>
constexpr Integer<T> operator*(Integer<T> lhs, S const& rhs)
{
    return lhs *= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator*(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) *= rhs;
}

template <typename T, typename S>
constexpr Integer<T> operator%(Integer<T> lhs, S const& rhs)
{
    return lhs %= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator%(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) %= rhs;
}

template <typename T, typename S>
constexpr Integer<T> operator<<(Integer<T> lhs, S const& rhs)
{
    return lhs <<= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator<<(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) <<= rhs;
}

template <typename T, typename S>
constexpr Integer<T> operator>>(Integer<T> lhs, S const& rhs)
{
    return lhs >>= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator>>(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) >>= rhs;
}

template <typename T, typename S>
constexpr Integer<T> operator&(Integer<T> lhs, S const& rhs)
{
    return lhs &= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator&(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) &= rhs;
}

template <typename T, typename S>
constexpr Integer<T> operator|(Integer<T> lhs, S const& rhs)
{
    return lhs |= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator|(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) |= rhs;
}

template <typename T, typename S>
constexpr Integer<T> operator^(Integer<T> lhs, S const& rhs)
{
    return lhs ^= Integer<T>(rhs);
}

template <typename T, typename S>
constexpr Integer<T> operator^(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>(lhs) ^= rhs;
}

template <typename T>
inline constexpr std::ostream& operator<<(std::ostream& os, Integer<T> const& n)
{
    return os << n.toString();
}

template <typename T>
constexpr Integer<T> gcd(Integer<T> const& a, Integer<T> const& b)
{
    if (a.isNan() || b.isNan() || a.isInfinity() || b.isInfinity())
        return Integer<T>::nan();

    if (a < 0)
        return gcd(a.abs(), b);

    if (b < 0)
        return gcd(a, b.abs());

    if (a < b)
        return gcd(b, a);

    if (!b)
        return a;

    return gcd(b, a % b);
}

template <typename T, typename S>
constexpr Integer<T> gcd(Integer<T> const& a, S const& b)
{
    return gcd(a, Integer<T>(b));
}

template <typename T, typename S>
constexpr Integer<T> gcd(S const& a, Integer<T> const& b)
{
    return gcd(Integer<T>(a), b);
}

template <typename T>
constexpr Integer<T> factorial(Integer<T> const& n)
{
    if (n.isNan() || n.isInfinity())
        return n;

    assert(n >= 0);

    if (n == 0)
        return Integer<T>(1);

    return n * factorial(n - 1);
}

template <typename T>
constexpr Integer<T> pow(Integer<T> base, Integer<T> exp)
{
    assert(exp >= 0);

    if (base.isInfinity() || base.isNan())
        return base;

    if (exp.isNan())
        return exp;

    if (exp.isInfinity())
        return exp;

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
constexpr Integer<T> pow(Integer<T> const& base, S const& exp)
{
    return pow(base, Integer<T>(exp));
}

template <typename T, typename S>
constexpr Integer<T> pow(S const& base, Integer<T> const& exp)
{
    return pow(Integer<T>(base), exp);
}

template <typename T>
constexpr Integer<T> powm(Integer<T> base, Integer<T> exp, Integer<T> const& mod)
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
constexpr Integer<T> powm(Integer<T> const& base, S const& exp, U const& mod)
{
    return powm(base, Integer<T>(exp), Integer<T>(mod));
}

template <typename T>
constexpr Integer<T> abs(Integer<T> const& n)
{
    return n.abs();
}

template <typename T>
constexpr std::pair<Integer<T>, Integer<T> > computeQr(Integer<T> const& dividend, Integer<T> const& divisor)
{
    if (!divisor)
        return {Integer<T>::nan(), Integer<T>::nan()};
    else if (!dividend)
        return {Integer<T>{0}, Integer<T>{0}};
    else if (divisor.abs() > dividend.abs())
        return {Integer<T>{0}, dividend};
    else if (dividend < 0 && divisor < 0)
    {
        auto qr{computeQr(-dividend, -divisor)};
        ++qr.first;
        qr.second += divisor;

        return qr;
    }
    else if (dividend > 0 && divisor < 0)
    {
        auto qr{computeQr(dividend, -divisor)};

        qr.first = -qr.first - 1;
        qr.second += divisor;

        return qr;
    }
    else if (dividend < 0 && divisor > 0)
    {
        auto qr{computeQr(-dividend, divisor)};

        qr.first = -qr.first - 1;
        qr.second = divisor - qr.second;

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
constexpr std::pair<Integer<T>, Integer<T> > computeQr(Integer<T> const& dividend, S const& divisor)
{
    return computeQr(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
constexpr std::pair<Integer<T>, Integer<T> > computeQr(S const& dividend, Integer<T> const& divisor)
{
    return computeQr(Integer<T>(dividend), divisor);
}

template <typename T>
constexpr Integer<T> computeQuotient(Integer<T> const& dividend, Integer<T> const& divisor)
{
    if (!divisor)
        return Integer<T>::nan();
    else if (!dividend)
        return Integer<T>{0};
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
constexpr std::pair<Integer<T>, Integer<T> > computeQuotient(Integer<T> const& dividend, S const& divisor)
{
    return computeQuotient(dividend, Integer<T>(divisor));
}

template <typename T, typename S>
constexpr std::pair<Integer<T>, Integer<T> > computeQuotient(S const& dividend, Integer<T> const& divisor)
{
    return computeQuotient(Integer<T>(dividend), divisor);
}

using Integerc = Integer<unsigned char>;
using Integers = Integer<unsigned short>;
using Integeri = Integer<unsigned int>;
using Integerl = Integer<unsigned long>;
using Integerll= Integer<unsigned long long>;

#endif // INTEGER_H
