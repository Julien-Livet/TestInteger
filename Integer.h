#ifndef INTEGER
#define INTEGER

/*!
 *  \file Integer.h
 *  \brief Provide a class to manage large integer numbers
 *  \author Julien LIVET
 *  \version 1.0
 *  \date 28/12/2024
 */

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <random>
#include <sstream>
#include <vector>

#if __cplusplus >= 202002L
#include <format>
#endif

#define LONGEST_TYPE (unsigned long long)

template <typename T, typename Enable = void>
class Integer;

template <typename T>
class Integer<T, typename std::enable_if<std::is_unsigned<T>::value>::type>
{
    public:
        Integer() : bits_{0}
        {
        }

        template <typename S>
        explicit Integer(S n) : isPositive_{n >= 0}
        {
            bits_.reserve(std::max(LONGEST_TYPE{1}, sizeof(S) / sizeof(T)));

            if (n < 0)
                n = -n;

            if (sizeof(T) == sizeof(S))
                bits_.emplace_back(n);
            else
            {
                auto const shift{LONGEST_TYPE{1} << std::min(sizeof(T), sizeof(S)) * 8};

                for (size_t i{0}; i < bits_.capacity(); ++i)
                {
                    bits_.emplace_back(n % shift);
                    n /= shift;
                }

                std::reverse(bits_.begin(), bits_.end());
            }
        }

        explicit Integer(std::vector<T> const& bits, bool isPositive = true) : isPositive_{isPositive}, bits_{bits}
        {
            if (bits_.empty())
                bits_.emplace_back(T{0});
        }

        explicit Integer(char const* n, size_t base = 10) : Integer{std::string{n}, base}
        {
        }

        explicit Integer(std::string n, size_t base = 10)
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

            std::string const str{it, n.end()};

            if (str == "nan")
                setNan();
            else if (str == "inf")
                setInfinity();
            else
            {
                auto const isPositive{isPositive_};

                if (base == 2)
                {
                    if (*it == 'b')
                        ++it;
                    else if (n.substr(0, 2) == "0b")
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
                    if (*it == 'o')
                        ++it;
                    else if (n.substr(0, 2) == "0o")
                        it += 2;

                    auto otherIt{n.rbegin()};
                    Integer exp{0};

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '7')
                        {
                            *this += (*otherIt - '0') * pow(base, exp);
                            ++exp;
                        }

                        ++otherIt;
                    }
                }
                else if (base <= 10)
                {
                    auto otherIt{n.rbegin()};
                    Integer exp{0};

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= static_cast<char>('0' + base))
                        {
                            *this += (*otherIt - '0') * pow(base, exp);
                            ++exp;
                        }

                        ++otherIt;
                    }
                }
                else if (base < 16)
                {
                    auto otherIt{n.rbegin()};
                    Integer exp{0};

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '9')
                        {
                            *this += (*otherIt - '0') * pow(base, exp);
                            ++exp;
                        }
                        else if ('a' <= std::tolower(*otherIt) && std::tolower(*otherIt) <= static_cast<char>('a' + base - 10))
                        {
                            *this += (*otherIt - 'a' + 10) * pow(base, exp);
                            ++exp;
                        }

                        ++otherIt;
                    }
                }
                else if (base == 16)
                {
                    if (*it == 'x')
                        ++it;
                    else if (n.substr(0, 2) == "0x")
                        it += 2;

                    auto otherIt{n.rbegin()};
                    Integer exp{0};

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '9')
                        {
                            *this += (*otherIt - '0') * pow(base, exp);
                            ++exp;
                        }
                        else if ('a' <= std::tolower(*otherIt) && std::tolower(*otherIt) <= 'f')
                        {
                            *this += (*otherIt - 'a' + 10) * pow(base, exp);
                            ++exp;
                        }

                        ++otherIt;
                    }
                }
                else// if (base <= 62)
                {
                    auto otherIt{n.rbegin()};
                    Integer exp{0};

                    while (otherIt.base() != it)
                    {
                        if ('0' <= *otherIt && *otherIt <= '9')
                        {
                            *this += (*otherIt - '0') * pow(base, exp);
                            ++exp;
                        }
                        else if ('a' <= *otherIt && *otherIt <= 'z')
                        {
                            *this += (*otherIt - 'a' + 10) * pow(base, exp);
                            ++exp;
                        }
                        else if ('A' <= *otherIt && *otherIt <= 'Z')
                        {
                            *this += (*otherIt - 'A' + 36) * pow(base, exp);
                            ++exp;
                        }

                        ++otherIt;
                    }
                }

                isPositive_ = isPositive;
            }
        }

        bool isPositive() const noexcept
        {
            return isPositive_;
        }

        bool isNegative() const noexcept
        {
            return !isPositive_;
        }

        auto const& bits() const noexcept
        {
            return bits_;
        }

        void invert() noexcept
        {
            for (size_t i{0}; i < bits_.size(); ++i)
                bits_[i] = ~bits_[i];
        }

        Integer& operator*=(Integer const& other)
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
                {/*
                    Integer result{0};

                    while (other > 0)
                    {
                        if (other & 1)
                            result += *this;

                        *this <<= 1;
                        other >>= 1;
                    }

                    *this = result;*/

                    if (*this == bits_.back() && other == other.bits_.back())
                    {
                        auto const n{LONGEST_TYPE{bits_.back()} * LONGEST_TYPE{other.bits_.back()}};

                        if (n / LONGEST_TYPE{other.bits_.back()} == LONGEST_TYPE{bits_.back()})
                            *this = n;
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
                            size_t const m{n / 2};
                            Integer x0, x1, y0, y1;
                            x0.bits_ = std::vector<T>(1, T{0});
                            x0.bits_.back() = ~(T{1} >> (sizeof(T) * 8 - m)) & bits_.back();
                            x1.bits_ = std::vector<T>(1, T{0});
                            x1.bits_.back() = ((~(T{1} >> (sizeof(T) * 8 - m)) << m) >> m) & bits_.back();
                            y0.bits_ = std::vector<T>(1, T{0});
                            y0.bits_.back() = ~(T{1} >> (sizeof(T) * 8 - m)) & other.bits_.back();
                            y1.bits_ = std::vector<T>(1, T{0});
                            y1.bits_.back() = ((~(T{1} >> (sizeof(T) * 8 - m)) << m) >> m) & other.bits_.back();

                            auto const z0{x0 * y0};
                            auto const z1{x1 * y0 + x0 * y1};
                            auto const z2{x1 * y1};

                            //xy = z2 * 2^(2 * m) + z1 * 2^m + z0
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
                        //Karatsuba algorithm
                        //x = x1 * 2^m + x0
                        //y = y1 * 2^m + y0
                        size_t n{std::max(number(), other.number()) / (sizeof(T) * 8)};
                        if (std::max(number(), other.number()) % (sizeof(T) * 8))
                            ++n;
                        size_t const m{n / 2};
                        std::cout << "*this ";
                        for (auto const& b : bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        Integer x0, x1, y0, y1;
                        x0.bits_ = std::vector<T>(m, T{0});
                        std::copy(bits_.rbegin(), bits_.rbegin() + m, x0.bits_.rbegin());
                        x1.bits_ = std::vector<T>(m, T{0});
                        std::copy(bits_.rbegin() + m, bits_.rend(), x1.bits_.rbegin());
                        std::cout << "x1 x0 ";
                        for (auto const& b : x1.bits_)
                            std::cout << (unsigned long long)b << " ";
                        for (auto const& b : x0.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << "other ";
                        for (auto const& b : other.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        y0.bits_ = std::vector<T>(m, T{0});
                        std::copy(other.bits_.rbegin(), other.bits_.rbegin() + m, y0.bits_.rbegin());
                        y1.bits_ = std::vector<T>(m, T{0});
                        std::copy(other.bits_.rbegin() + m, other.bits_.rend(), y1.bits_.rbegin());
                        std::cout << "y1 y0 ";
                        for (auto const& b : y1.bits_)
                            std::cout << (unsigned long long)b << " ";
                        for (auto const& b : y0.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;

                        auto const z0{x0 * y0};
                        auto const z1{x1 * y0 + x0 * y1};
                        auto const z2{x1 * y1};

                        //xy = z2 * 2^(2 * m) + z1 * 2^m + z0
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

            return *this;
        }

        Integer& operator+=(Integer const& other)
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
                        T const bit_a{(i < a.size()) ? a[a.size() - 1 - i] : T{0}};
                        T const bit_b{(i < b.size()) ? b[b.size() - 1 - i] : T{0}};

                        auto const sum{static_cast<T>(bit_a + bit_b + carry)};
                        carry = (sum - bit_b < bit_a + carry);

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
                        T const bit_a{(i < a.size()) ? a[a.size() - 1 - i] : T{0}};
                        T const bit_b{(i < b.size()) ? b[b.size() - 1 - i] : T{0}};

                        T const bit_result{bit_a - bit_b - borrow};
                        borrow = (bit_result + bit_b > bit_a - borrow);

                        result.emplace_back(bit_result);
                    }

                    std::reverse(result.begin(), result.end());

                    bits_ = result;
                }
            }

            return *this;
        }

        Integer& operator-=(Integer const& other)
        {
            return *this += -other;
        }

        Integer& operator/=(Integer other)
        {
            if (other.isNegative())
            {
                *this = -*this;

                return *this /= -other;
            }

            auto const n{*this};

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
                    if (~T{0} >= other.abs() && bits_.size() == 1)
                        *this = Integer{bits_.back() / other.template cast<unsigned long long>()};
                    else
                    {/**
                        auto temp_d{other.abs()};
                        Integer temp_q{1};
                        Integer q{0};
                        auto r{abs()};

                        while (temp_d <= r)
                        {
                            temp_d <<= 1;
                            temp_q <<= 1;
                        }

                        while (temp_q > 0)
                        {
                            if (r >= temp_d)
                            {
                                r -= temp_d;
                                q |= temp_q;
                            }

                            temp_d >>= 1;
                            temp_q >>= 1;
                        }

                        *this = q;**//*
                        auto const d{other.abs()};
                        Integer q{0};
                        auto r{abs()};
                        Integer p{1};

                        while (other <= r && p <= r)
                        {
                            other <<= 1;
                            p <<= 1;
                        }

                        while (p != 0)
                        {
                            other >>= 1;
                            p >>= 1;

                            if (r >= other)
                            {
                                r -= other;
                                q += p;
                            }
                        }
                        std::cout << "q ";
                        for (auto const& b : q.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        auto const r_{*this - q * d};
                        std::cout << "r_ ";
                        for (auto const& b : r_.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        *this = q;*/

                        Integer q{0};
                        Integer r{0};
                        size_t const iMax{number() - 1};

                        for (size_t i{iMax}; i <= iMax; --i)
                        {
                            r <<= 1;
                            r.setBit(0, bit(i));

                            if (r >= other.abs())
                            {
                                r -= other.abs();
                                q.setBit(i, true);
                            }
                        }

                        std::cout << abs().toString(2) << std::endl;
                        std::cout << (q * other.abs() + r).toString(2) << std::endl;
                        //assert(abs() == q * other.abs() + r);

                        *this = q;
                    }
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

        Integer& operator%=(Integer const& other)
        {
            if (!other || other.isNan() || other.isInfinity())
                setNan();
            else
            {
                if ((isPositive_ && other.isPositive_) ||
                    (!isPositive_ && !other.isPositive_))
                {
                    if (~T{0} >= other.abs() && bits_.size() == 1)
                    {
                        auto const isPositive{isPositive_};

                        *this = Integer{bits_.back() % other.abs().template cast<T>()};

                        isPositive_ = isPositive;
                    }
                    else
                    {/*
                        auto const n{*this / other};
                        std::cout << "this ";
                        for (auto const& b : bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << "other ";
                        for (auto const& b : other.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << "n ";
                        for (auto const& b : n.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        auto const n1{other * n};
                        std::cout << "n1 ";
                        for (auto const& b : n1.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        auto const n2{*this - n1};
                        std::cout << "n2 ";
                        for (auto const& b : n2.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        *this -= other * n;*//*
                        std::cout << "q ";
                        for (auto const& b : (*this / other).bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << "r ";
                        for (auto const& b : (*this - other * (*this / other)).bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        *this -= other * (*this / other);*/
/**
                        auto temp_d{other.abs()};
                        Integer temp_q{1};
                        Integer q{0};
                        auto r{abs()};

                        while (temp_d <= r)
                        {
                            temp_d <<= 1;
                            temp_q <<= 1;
                        }

                        while (temp_q > 0)
                        {
                            if (r >= temp_d)
                            {
                                r -= temp_d;
                                q |= temp_q;
                            }

                            temp_d >>= 1;
                            temp_q >>= 1;
                        }
                        std::cout << "this ";
                        for (auto const& b : bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << "other ";
                        for (auto const& b : other.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << "q ";
                        for (auto const& b : q.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        auto const num1{q * other};
                        std::cout << "num1 ";
                        for (auto const& b : num1.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        auto const num2{*this - num1};
                        std::cout << "num2 ";
                        for (auto const& b : num2.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << "r ";
                        for (auto const& b : r.bits_)
                            std::cout << (unsigned long long)b << " ";
                        std::cout << std::endl;
                        std::cout << (temp_d <= r) << std::endl;

                        *this = r;**/
                        Integer q{0};
                        Integer r{0};
                        size_t const iMax{number() - 1};

                        for (size_t i{iMax}; i <= iMax; --i)
                        {
                            r <<= 1;
                            r.setBit(0, bit(i));

                            if (r >= other.abs())
                            {
                                r -= other.abs();
                                q.setBit(i, true);
                            }
                        }

                        //assert(abs() == q * other.abs() + r);

                        *this = r;
                    }
                }
                else if (isPositive_ && !other.isPositive_)
                {
                    *this %= -other;
                    *this -= -other;
                }
                else if (!isPositive_ && other.isPositive_)
                {
                    *this %= -other;
                    *this += other;
                }
            }

            if (!(abs() < other.abs())){
                assert(abs() < other.abs());
            }

            return *this;
        }

        Integer& operator<<=(Integer other)
        {
            assert(other >= 0);

            if (!*this)
                return *this;

            auto const s{static_cast<unsigned short>(sizeof(T) * 8)};
            auto const n{other / s};

            std::vector<T> const v(n.template cast<unsigned long long>(), T{0});

            bits_.insert(bits_.end(), v.begin(), v.end());

            other -= n * s;

            std::vector<T> bits(bits_.size() + 1, T{0});

            std::copy(bits_.rbegin(), bits_.rend(), bits.rbegin());

            bits_ = bits;

            auto const shift{other.template cast<unsigned long long>()};

            for (auto it{bits_.begin() + 1}; it != bits_.end(); ++it)
            {
                if ((*it >> (sizeof(T) * 8 - shift)))
                    *(it - 1) |= (*it >> (sizeof(T) * 8 - shift));

                *it <<= shift;
            }

            if (!bits_.front())
            {
                std::copy(bits_.begin() + 1, bits_.end(), bits_.begin());
                bits_.pop_back();
            }

            return *this;
        }

        Integer& operator>>=(Integer other)
        {
            assert(other >= 0);

            if (!*this)
                return *this;

            auto const s{static_cast<unsigned short>(sizeof(T) * 8)};
            auto const n{other / s};

            bits_.resize(bits_.size() - n.template cast<unsigned long long>());

            other -= n * s;

            auto const shift{other.template cast<unsigned long long>()};

            for (auto it{bits_.rbegin()}; it != bits_.rend(); ++it)
            {
                *it >>= shift;

                if (it != bits_.rend() - 1 && (*(it + 1) & ((1 << shift) - 1)))
                    *it |= (*(it +  1) & ((1 << shift) - 1)) << (sizeof(T) * 8 - shift);
            }

            if (bits_.empty())
                bits_.emplace_back(T{0});

            return *this;
        }

        bool operator>=(Integer const& other) const
        {
            return operator>(other) || operator==(other);
        }

        bool operator>(Integer const& other) const
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

        bool operator<=(Integer const& other) const
        {
            return operator<(other) || operator==(other);
        }

        bool operator<(Integer const& other) const
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

        bool operator==(Integer const& other) const noexcept
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
        bool operator==(S const& other) const
        {
            return *this == Integer{other};
        }

        bool operator!=(Integer const& other) const
        {
            return !operator==(other);
        }

        template <typename S>
        bool operator!=(S const& other) const
        {
            return *this != Integer{other};
        }

        Integer operator-() const
        {
            auto x{*this};

            x.isPositive_ = !x.isPositive_;

            return x;
        }

        Integer operator~() const
        {
            auto x{*this};

            x.invert();

            return x;
        }

        operator bool() const noexcept
        {
            return !!*this;
        }

        bool operator!() const noexcept
        {
            for (auto const& b : bits_)
            {
                if (b)
                    return false;
            }

            return true;
        }

        Integer& operator--()
        {
            return *this -= 1;
        }

        Integer operator--(int)
        {
            auto x{*this};

            operator--();

            return x;
        }

        Integer& operator++()
        {
            return *this += 1;
        }

        Integer operator++(int)
        {
            auto x{*this};

            operator++();

            return x;
        }

        template <typename S>
        Integer& operator+=(S const& other)
        {
            return *this += Integer{other};
        }

        template <typename S>
        Integer& operator-=(S const& other)
        {
            return *this -= Integer{other};
        }

        template <typename S>
        Integer& operator/=(S const& other)
        {
            return *this /= Integer{other};
        }

        template <typename S>
        Integer& operator*=(S const& other)
        {
            return *this *= Integer{other};
        }

        template <typename S>
        Integer& operator%=(S const& other)
        {
            return *this %= Integer{other};
        }

        template <typename S>
        Integer& operator>>=(S const& other)
        {
            return *this >>= Integer{other};
        }

        template <typename S>
        Integer& operator<<=(S const& other)
        {
            return *this <<= Integer{other};
        }

        Integer& operator&=(Integer const& other)
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

        Integer& operator|=(Integer const& other)
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
        Integer& operator|=(S const& other)
        {
            return *this |= Integer{other};
        }

        template <typename S>
        Integer& operator&=(S const& other)
        {
            return *this &= Integer{other};
        }

        template <typename S>
        Integer& operator^=(S const& other)
        {
            return *this ^= Integer{other};
        }

        Integer& operator^=(Integer const& other)
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
        Integer& operator=(S const& other)
        {
            return *this = Integer{other};
        }

        std::string toString(size_t base = 10) const
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
                    auto number{abs()};

                    if (!number)
                        s = "0";

                    while (number)
                    {
                        auto const tmp{number % 8};
                        s = std::to_string(tmp.template cast<short>()) + s;
                        number /= 8;
                    }
                }
#endif

                s = "0o" + s;
            }
            else if (base == 10)
            {
                auto number{abs()};

                if (bits_.size() == 1)
                    s = std::to_string(bits_.back());
                else
                {
                    if (!number)
                        s = "0";

                    auto const n{1};//auto const n{static_cast<T>(std::log10(~T{0}))};
                    Integer const b{pow(Integer{10}, n)};
                    //std::cout << "ho " << (unsigned long long)n << " " << (b < ~T{0}) << std::endl;

                    while (number)
                    {
                        std::cout << "number ";
                        for (auto const& bit : number.bits_)
                            std::cout << (unsigned long long)bit << " ";
                        std::cout << std::endl;
                        std::cout << "b ";
                        for (auto const& bit : b.bits_)
                            std::cout << (unsigned long long)bit << " ";
                        std::cout << std::endl;
                        auto const tmp{number % b};
                        std::cout << "tmp ";
                        for (auto const& bit : tmp.bits_)
                            std::cout << (unsigned long long)bit << " ";
                        std::cout << std::endl;
                        //std::cout << "he " << (unsigned long long)tmp.template cast<T>() << std::endl;
                        //std::cout << (unsigned long long)number.bits_.back() << " " << std::endl;
                        std::ostringstream oss;
                        oss << std::setw(n) << std::setfill('0') << tmp.template cast<unsigned long long>();
                        s = oss.str() + s;
                        //std::cout << "a " << " " << (number < b) << std::endl;
                        number /= b;
                        std::cout << "number ";
                        for (auto const& bit : number.bits_)
                            std::cout << (unsigned long long)bit << " ";
                        std::cout << std::endl;
                        //std::cout << (unsigned long long)number.bits_.back() << " " << std::endl;
                        //std::string s;std::getline(std::cin, s);
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
                auto number{abs()};

                if (bits_.size() == 1)
                    s = std::to_string(bits_.back());
                else
                {
                    if (!number)
                        s = "0";

                    while (number)
                    {
                        auto tmp{number};
                        tmp %= static_cast<unsigned char>(base);
                        s = std::to_string(tmp.template cast<short>()) + s;
                        number /= static_cast<unsigned char>(base);
                    }
                }
            }
            else if (base == 16)
            {
                auto number{abs()};
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
                        auto const tmp{number % 16};
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
                auto number{abs()};

                if (!number)
                    s = "0";

                while (number)
                {
                    auto const tmp{number % 62};
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

        operator char() const noexcept
        {
            return cast<char>();
        }

        operator unsigned char() const noexcept
        {
            return cast<unsigned char>();
        }

        operator short() const noexcept
        {
            return cast<short>();
        }

        template <typename S>
        S cast() const noexcept
        {
            S n{0};

            size_t const iMax{std::min(std::max(LONGEST_TYPE{1}, sizeof(S) / sizeof(T)), bits_.size())};
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

        operator unsigned short() const noexcept
        {
            return cast<unsigned short>();
        }

        operator int() const noexcept
        {
            return cast<int>();
        }

        operator unsigned int() const noexcept
        {
            return cast<unsigned int>();
        }

        operator long() const noexcept
        {
            return cast<long>();
        }

        operator unsigned long() const noexcept
        {
            return cast<unsigned long>();
        }

        operator long long() const noexcept
        {
            return cast<long long>();
        }

        operator unsigned long long() const noexcept
        {
            return cast<unsigned long long>();
        }

        bool isNan() const noexcept
        {
            return isNan_;
        }

        void setNan() noexcept
        {
            isNan_ = true;
            isInfinity_ = false;
            bits_.clear();
        }

        bool isInfinity() const noexcept
        {
            return isInfinity_;
        }

        void setInfinity() noexcept
        {
            isNan_ = false;
            isInfinity_ = true;
            bits_.clear();
        }

        Integer abs() const
        {
            if (isNegative())
                return -*this;

            return *this;
        }

        size_t precision() const noexcept
        {
            return bits_.size();
        }

        void setPrecision(size_t precision)
        {
            assert(precision);

            std::vector<T> bits(precision, T{0});

            std::copy(bits_.rbegin(), bits_.rbegin() + std::min(bits_.size(), precision), bits.rbegin());

            bits_ = bits;
        }

        template <typename URNG>
        void setRandom(URNG& g) noexcept
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

        bool isPrime(size_t reps = 50) const
        {
            if (*this < 2)
                return false;

            if (*this != 2 && *this % 2 == 0)
                return false;

            auto s{*this - 1};

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
            auto const number{*this - 1};

            for (size_t i{0}; i < reps; ++i)
            {
                auto a{*this};
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

        void setPositive()
        {
            isPositive_ = true;
        }

        void setNegative()
        {
            isPositive_ = false;
        }

        static Integer nan()
        {
            Integer n;
            n.setNan();

            return n;
        }

        static Integer infinity()
        {
            static Integer n;
            n.setInfinity();

            return n;
        }

        bool bit(size_t n) const noexcept
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

        void setBit(size_t n, bool bit)
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

        size_t count() const noexcept
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

        size_t number() const noexcept
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

        bool isEven() const noexcept
        {
            if (bits_.empty())
                return false;

            return !(bits_.back() & 1);
        }

        bool isOdd() const noexcept
        {
            if (bits_.empty())
                return false;

            return bits_.back() & 1;
        }

        template <typename S>
        bool fits() const
        {
            return (*this == this->template cast<S>());
        }

        Integer sign() const
        {
            if (*this < 0)
                return Integer{-1};

            return Integer{1};
        }

        void setSign(Integer const& other) noexcept
        {
            isPositive_ = other.isPositive_;
        }

        Integer previousPrime() const
        {
            if (isNan())
                return *this;

            if (isInfinity() || *this < 2)
                return nan();

            if (*this == 2)
                return 2;
            else if (*this == 3)
                return 2;

            Integer n{*this - 2};

            while (!n.isPrime())
                n -= 2;
        }

        Integer nextPrime() const
        {
            if (isNan())
                return *this;

            if (*this < 2)
                return 2;
            else if (*this == 2)
                return 3;
            else if (isInfinity())
                return nan();

            Integer n{*this + 2};

            while (!n.isPrime())
                n += 2;
        }

        size_t size() const noexcept
        {
            return bits_.size();
        }

    private:
        bool isPositive_{true};
        std::vector<T> bits_;
        bool isNan_{false};
        bool isInfinity_{false};
};

template <typename T>
Integer<T> operator*(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs *= rhs;
}

template <typename T>
Integer<T> operator+(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs += rhs;
}

template <typename T>
Integer<T> operator-(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs -= rhs;
}

template <typename T>
Integer<T> operator/(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs /= rhs;
}

template <typename T>
Integer<T> operator%(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs %= rhs;
}

template <typename T>
Integer<T> operator&(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs &= rhs;
}

template <typename T>
Integer<T> operator|(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs |= rhs;
}

template <typename T>
Integer<T> operator^(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs ^= rhs;
}

template <typename T>
Integer<T> operator<<(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs <<= rhs;
}

template <typename T>
Integer<T> operator>>(Integer<T> lhs, Integer<T> const& rhs)
{
    return lhs >>= rhs;
}

template <typename T, typename S>
bool operator>(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator>(Integer<T>{rhs});
}

template <typename T, typename S>
bool operator>(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator<(Integer<T>{lhs});
}

template <typename T, typename S>
bool operator>=(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator>=(Integer<T>{rhs});
}

template <typename T, typename S>
bool operator>=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator<=(Integer<T>{lhs});
}

template <typename T, typename S>
bool operator<(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator<(Integer<T>{rhs});
}

template <typename T, typename S>
bool operator<(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator>(Integer<T>{lhs});
}

template <typename T, typename S>
bool operator<=(Integer<T> const& lhs, S const& rhs) noexcept
{
    return lhs.operator<=(Integer<T>{rhs});
}

template <typename T, typename S>
bool operator<=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator>=(Integer<T>{lhs});
}

template <typename T, typename S>
bool operator==(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator==(Integer<T>{lhs});
}

template <typename T, typename S>
bool operator!=(S const& lhs, Integer<T> const& rhs) noexcept
{
    return rhs.operator!=(Integer<T>{lhs});
}

template <typename T, typename S>
Integer<T> operator+(Integer<T> lhs, S const& rhs)
{
    return lhs += Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator+(S const& lhs, Integer<T> rhs)
{
    return rhs += Integer<T>{lhs};
}

template <typename T, typename S>
Integer<T> operator-(Integer<T> lhs, S const& rhs)
{
    return lhs -= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator-(S const& lhs, Integer<T> rhs)
{
    return -(rhs -= Integer<T>{lhs});
}

template <typename T, typename S>
Integer<T> operator/(Integer<T> lhs, S const& rhs)
{
    return lhs /= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator/(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} /= rhs;
}

template <typename T, typename S>
Integer<T> operator*(Integer<T> lhs, S const& rhs)
{
    return lhs *= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator*(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} *= rhs;
}

template <typename T, typename S>
Integer<T> operator%(Integer<T> lhs, S const& rhs)
{
    return lhs %= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator%(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} %= rhs;
}

template <typename T, typename S>
Integer<T> operator<<(Integer<T> lhs, S const& rhs)
{
    return lhs <<= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator<<(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} <<= rhs;
}

template <typename T, typename S>
Integer<T> operator>>(Integer<T> lhs, S const& rhs)
{
    return lhs >>= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator>>(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} >>= rhs;
}

template <typename T, typename S>
Integer<T> operator&(Integer<T> lhs, S const& rhs)
{
    return lhs &= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator&(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} &= rhs;
}

template <typename T, typename S>
Integer<T> operator|(Integer<T> lhs, S const& rhs)
{
    return lhs |= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator|(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} |= rhs;
}

template <typename T, typename S>
Integer<T> operator^(Integer<T> lhs, S const& rhs)
{
    return lhs ^= Integer<T>{rhs};
}

template <typename T, typename S>
Integer<T> operator^(S const& lhs, Integer<T> const& rhs)
{
    return Integer<T>{lhs} ^= rhs;
}

template <typename T>
inline std::ostream& operator<<(std::ostream& os, Integer<T> const& n)
{
    return os << n.toString();
}

template <typename T>
Integer<T> gcd(Integer<T> const& a, Integer<T> const& b)
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
Integer<T> gcd(Integer<T> const& a, S const& b)
{
    return gcd(a, Integer<T>{b});
}

template <typename T, typename S>
Integer<T> gcd(S const& a, Integer<T> const& b)
{
    return gcd(Integer<T>{a}, b);
}

template <typename T>
Integer<T> factorial(Integer<T> const& n)
{
    if (n.isNan() || n.isInfinity())
        return n;

    assert(n >= 0);

    if (n == 0)
        return Integer<T>{1};

    return n * factorial(n - 1);
}

template <typename T>
Integer<T> pow(Integer<T> base, Integer<T> exp)
{
    assert(exp >= 0);

    if (base.isInfinity() || base.isNan())
        return base;

    if (exp.isNan())
        return exp;

    if (exp.isInfinity())
        return exp;

    if (base == 2)
        return Integer<T>{1} << exp;

    Integer<T> result{1};

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
Integer<T> pow(Integer<T> const& base, S const& exp)
{
    return pow(base, Integer<T>{exp});
}

template <typename T, typename S>
Integer<T> pow(S const& base, Integer<T> const& exp)
{
    return pow(Integer<T>{base}, exp);
}

template <typename T>
Integer<T> powm(Integer<T> base, Integer<T> exp, Integer<T> const& mod)
{
    assert(exp >= 0);

    Integer<T> result{1};

    auto base_mod{base % mod};

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
Integer<T> powm(Integer<T> const& base, S const& exp, U const& mod)
{
    return powm(base, Integer<T>{exp}, Integer<T>{mod});
}

template <typename T>
Integer<T> abs(Integer<T> const& n)
{
    return n.abs();
}

using Integerc = Integer<unsigned char>;
using Integers = Integer<unsigned short>;
using Integeri = Integer<unsigned int>;
using Integerl = Integer<unsigned long>;
using Integerll= Integer<unsigned long long>;

#endif // INTEGER
