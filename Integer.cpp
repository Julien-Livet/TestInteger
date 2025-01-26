#include <atomic>
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

Integer::Integer(std::initializer_list<uintmax_t> const& bits, bool isPositive) : isPositive_{isPositive}, bits_{bits}
{
    if (bits_.empty())
        bits_.emplace_back(0);

    adjust();
}

#ifdef WITH_GMP
Integer::Integer(mpz_class const& n) : Integer(n.get_str(2), 2)
{
}
#endif

Integer::Integer(char const* n, size_t base) : Integer(std::string{n}, base)
{
}

Integer::Integer(std::string n, size_t base)
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
            if (std::tolower(*it) == 'b')
                ++it;
            else if (str.substr(0, 2) == "0b")
                it += 2;

            *this = 0;

            while (it != n.end())
            {
                if (*it == '1')
                {
                    *this <<= 1;
                    *this |= 1;
                }
                else if (*it == '0')
                    *this <<= 1;

                ++it;
            }
        }
        else if (base == 8)
        {
            if (str[0] == 'o')
                ++it;
            else if (str.substr(0, 2) == "0o")
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
            if (str[0] == 'x')
                ++it;
            else if (str.substr(0, 2) == "0x")
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

CONSTEXPR std::vector<uintmax_t> const& Integer::bits() const noexcept
{
    return bits_;
}

void Integer::invert() noexcept
{
    auto threadFunc
        {
            [this] (size_t start, size_t end) -> void
            {
                for (size_t i{start}; i < end; ++i)
                    this->bits_[i] = ~this->bits_[i];
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

    if (autoAdjust_)
        adjust();
}

CONSTEXPR bool Integer::isNegative() const noexcept
{
    return !isPositive_;
}

CONSTEXPR bool Integer::isPositive() const noexcept
{
    return isPositive_;
}

bool Integer::operator>=(Integer const& other) const
{
    return operator>(other) || operator==(other);
}

bool Integer::operator>(Integer const& other) const
{
    if (!isPositive_ && other.isPositive_)
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

    std::vector<uintmax_t> a(std::max(bits_.size(), other.bits_.size()), 0);
    std::vector<uintmax_t> b{a};

    std::copy(bits_.rbegin(), bits_.rend(), a.rbegin());
    std::copy(other.bits_.rbegin(), other.bits_.rend(), b.rbegin());

    auto const great{a > b};

    return isPositive_ ? great : !great;
}

bool Integer::operator<=(Integer const& other) const
{
    return operator<(other) || operator==(other);
}

bool Integer::operator<(Integer const& other) const
{
    if (isPositive_ && !other.isPositive_)
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

    std::vector<uintmax_t> a(std::max(bits_.size(), other.bits_.size()), 0);
    std::vector<uintmax_t> b{a};

    std::copy(bits_.rbegin(), bits_.rend(), a.rbegin());
    std::copy(other.bits_.rbegin(), other.bits_.rend(), b.rbegin());

    auto const less{a < b};

    return isPositive_ ? less : !less;
}

CONSTEXPR bool Integer::operator==(Integer const& other) const noexcept
{
    if (isNan() && other.isNan())
        return true;
    else if (isNan() || other.isNan())
        return false;
    else if (isInfinity() && other.isInfinity())
        return (isPositive() == other.isPositive());
    else if (isInfinity() || other.isInfinity())
        return false;

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

CONSTEXPR bool Integer::operator!=(Integer const& other) const
{
    return !operator==(other);
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

bool Integer::operator!() const noexcept
{
    for (auto const& b : bits_)
    {
        if (b)
            return false;
    }

    return true;
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

std::string Integer::toString(size_t base) const
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
        s.reserve(s.size() + 2 + bits_.size() * sizeof(uintmax_t) * 8);

#if __cplusplus >= 202106L
        switch (sizeof(uintmax_t))
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

            for (size_t i{0}; i < sizeof(uintmax_t) * 8; ++i)
            {
                s = (b & 1 ? '1' : '0') + s;
                b >>= 1;
            }
        }
#endif

        s = "0b" + s;
    }
    else if (base == 8)
    {
#if __cplusplus >= 202106L
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
                s = std::to_string(tmp.cast<short>()) + s;
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

            auto const n{static_cast<uintmax_t>(std::log10(static_cast<uintmax_t>(~uintmax_t{0})))};
            uintmax_t const b(pow(uintmax_t{10}, n));

            while (number)
            {
                auto const tmp(number % b);
                std::ostringstream oss;
                oss << std::setw(n) << std::setfill('0') << tmp.cast<uintmax_t>();
                s = oss.str() + s;
                number /= b;
            }

            size_t i{0};

            while (s[i] == '0' && i != s.size())
                ++i;

            if (i == s.size())
                i = s.size() - 1;

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
                auto const tmp(number % static_cast<unsigned char>(base));
                s = std::to_string(tmp.cast<short>()) + s;
                number /= static_cast<unsigned char>(base);
            }
        }
    }
    else if (base == 16)
    {
        auto number(abs());
#if __cplusplus >= 202106L
        switch (sizeof(uintmax_t))
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
                    s = std::to_string(tmp.cast<short>()) + s;
                else
                    s = (char)('a' + tmp.cast<short>() - 10) + s;
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
            auto const tmp(number % base);
            if (number < 10)
                s = std::to_string(tmp.cast<short>()) + s;
            else if (number - 10 < 26)
                s = (char)('a' + tmp.cast<short>() - 10) + s;
            else
                s = (char)('A' + tmp.cast<short>() - 36) + s;
            number /= base;
        }
    }

    if (bits_.empty() && (!s.empty() && s.back() != '0'))
        s += '0';

    if (!isPositive_)
        s = '-' + s;

    return s;
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

size_t Integer::precision() const noexcept
{
    return bits_.size();
}

void Integer::setPrecision(size_t precision)
{
    assert(precision);

    if (precision == bits_.size())
        return;

    std::vector<uintmax_t> bits(precision, 0);

    std::copy(bits_.rbegin(), bits_.rbegin() + std::min(bits_.size(), precision), bits.rbegin());

    bits_ = bits;
}

int Integer::isPrime(size_t reps) const
{
    assert(reps);

    if (*this < 2)
        return 0;
    else if (*this == 2)
        return 2;
    else if(!(*this & 1))
        return 0;
    else if (this->template fits<unsigned int>())
    {
        auto p{std::equal_range(primes.begin(), primes.end(), this->cast<unsigned int>())};

        if (p.first != primes.end() && *p.first != this->cast<unsigned int>())
            --p.first;

        if (p.first != p.second && *p.first == this->cast<unsigned int>())
            return 2;
    }

    auto isPrimeDivisible
        {
            [](Integer const& n, unsigned int prime) -> bool
            {
                return !(n % prime);
            }
        };

    //Trial divisions

    {
        auto const sqrtLimit(sqrt(*this));
        std::atomic<bool> divisible(false);

        auto threadFunc
            {
                [this, &divisible, &sqrtLimit, &isPrimeDivisible] (size_t start, size_t end) -> void
                {
                    for (size_t i{start}; i < end && !divisible.load(); ++i)
                    {
                        if (primes[i] > sqrtLimit)
                            break;

                        if (isPrimeDivisible(*this, primes[i]))
                        {
                            divisible.store(true);
                            return;
                        }
                    }
                }
            };

        size_t const numThreads{std::thread::hardware_concurrency()};
        size_t const chunkSize{primes.size() / numThreads};
        std::vector<std::thread> threads;

        for (size_t i{0}; i < numThreads; ++i)
        {
            size_t const start{i * chunkSize};
            size_t const end{(i == numThreads - 1) ? primes.size() : (i + 1) * chunkSize};
            threads.emplace_back(threadFunc, start, end);
        }

        for (auto& t : threads)
            t.join();

        if (divisible.load())
            return 0;

        if (sqrtLimit < primes.back())
            return 2;
    }

    //Miller-Rabin tests

    Integer s(*this - 1);

    while (!(s & 1))
        s >>= 1;

    auto reduction
    {
        [] (Integer const& t, Integer const& R, Integer const& n, Integer const& n_) -> Integer
        {
            auto const m(((t % R) * n_) % R);
            auto const x((t + m * n) / R);

            if (x < n)
                return x;
            else
                return x - n;
        }
    };

    auto redmulmod
    {
        [&reduction] (Integer const& a, Integer b, Integer const& n,
                      Integer const& R, Integer const& n_, Integer const& R2modn) -> Integer
        {
            auto const reda(reduction(a * R2modn, R, n, n_));
            auto const redb(reduction(b * R2modn, R, n, n_));
            auto const redc(reduction(reda * redb, R, n, n_));

            return reduction(redc, R, n, n_);
        }
    };

    auto mulmod
    {
        [] (Integer const& a, Integer b, Integer const& m) -> bool
        {
            Integer x(0);
            Integer y(a % m);

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

    auto modulo
    {
        [&redmulmod] (Integer const& base, Integer e, Integer const& m,
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

    auto const number(*this - 1);
    auto const& m(*this);
    Integer R(Integer(1) << m.number());

    assert(R > m);

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

    auto const R2modm((R * R) % m);
    Integer R_, m_;
    auto const d(gcdExtended(R, -m, R_, m_));

    assert(d.abs() == 1);
    assert(R * R_ - m * m_ == d);

    if (d == -1)
    {
        R_ = -R_;
        m_ = -m_;
    }

    {
        std::atomic<bool> divisible(false);

        auto threadFunc
        {
            [this, &divisible, &modulo, &mulmod,
             &R, &m_, &R2modm, &number, &s]
            (size_t start, size_t end) -> void
            {
                for (size_t i{start}; i < end && !divisible.load(); ++i)
                {
                    auto a(*this);
                    a.template setRandom<std::random_device>();
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
                    {
                        divisible.store(true);
                        return;
                    }
                }
            }
        };

        size_t const numThreads{std::thread::hardware_concurrency()};
        size_t const chunkSize{reps / numThreads};
        std::vector<std::thread> threads;

        for (size_t i{0}; i < numThreads; ++i)
        {
            size_t const start{i * chunkSize};
            size_t const end{(i == numThreads - 1) ? reps : (i + 1) * chunkSize};
            threads.emplace_back(threadFunc, start, end);
        }

        for (auto& t : threads)
            t.join();

        if (divisible.load())
            return 0;
    }

    return 1;
}

CONSTEXPR void Integer::setPositive()
{
    isPositive_ = true;
}

CONSTEXPR void Integer::setNegative()
{
    isPositive_ = false;
}

bool Integer::bit(size_t n) const noexcept
{
    auto it{bits_.rbegin()};

    while (it != bits_.rend() && n > sizeof(uintmax_t) * 8)
    {
        n -= sizeof(uintmax_t) * 8;
        ++it;
    }

    if (it == bits_.rend())
        return false;

    return *it & (uintmax_t{1} << n);
}

void Integer::setBit(size_t n, bool bit)
{
    auto it{bits_.rbegin()};

    while (n > sizeof(uintmax_t) * 8)
    {
        n -= sizeof(uintmax_t) * 8;
        ++it;

        if (it == bits_.rend())
        {
            std::vector<uintmax_t> bits(bits_.size() + 1, 0);

            std::copy(bits_.rbegin(), bits_.rend(), bits.rbegin());

            bits_ = bits;

            it = bits_.rend() - 1;
        }
    }

    if (bit)
        *it |= uintmax_t{1} << n;
    else
        *it &= ~(uintmax_t{1} << n);
}

uintmax_t Integer::bits(size_t n) const noexcept
{
    if (n >= bits_.size())
        return 0;

    return bits_[bits_.size() - 1 - n];
}

void Integer::setBits(size_t n, uintmax_t const& bits)
{
    if (bits_.size() < n)
    {
        std::vector<uintmax_t> bits(n, 0);

        std::copy(bits_.rbegin(), bits_.rend(), bits.rbegin());

        bits_ = bits;
    }

    bits_[bits_.size() - 1 - n] = bits;

    if (autoAdjust_)
        adjust();
}

size_t Integer::count() const noexcept
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

bool Integer::isEven() const noexcept
{
    if (bits_.empty())
        return false;

    return !(bits_.back() & 1);
}

bool Integer::isOdd() const noexcept
{
    if (bits_.empty())
        return false;

    return bits_.back() & 1;
}

Integer Integer::sign() const
{
    if (*this < 0)
        return Integer(-1);

    return Integer(1);
}

CONSTEXPR void Integer::setSign(Integer const& other) noexcept
{
    isPositive_ = other.isPositive_;
}

Integer Integer::previousPrime() const
{
    if (isNan())
        return *this;
    else if (isInfinity() || *this < 2)
        return nan();
    else if (*this == 2)
        return Integer(2);
    else if (*this == 3)
        return Integer(2);
    else if (this->template fits<unsigned int>())
    {
        auto p{std::equal_range(primes.begin(), primes.end(), this->cast<unsigned int>())};

        if (p.first != primes.end() && *p.first != this->cast<unsigned int>())
            --p.first;

        if (p.first != p.second && p.second != primes.end())
        {
            if (*p.first == this->cast<unsigned int>())
                return Integer(*(p.first - 1));
            else
                return Integer(*p.first);
        }
    }

    auto n(*this);

    if (n % 2)
        n -= 2;
    else
        --n;

    while (!n.isPrime())
        n -= 2;

    return n;
}

Integer Integer::nextPrime() const
{
    if (isNan())
        return *this;
    else if (*this < 2)
        return Integer(2);
    else if (*this == 2)
        return Integer(3);
    else if (isInfinity())
        return nan();
    else if (this->template fits<unsigned int>())
    {
        auto p{std::equal_range(primes.begin(), primes.end(), this->cast<unsigned int>())};

        if (p.first != primes.end() && *p.first != this->cast<unsigned int>())
            --p.first;

        if (p.first != p.second && p.second != primes.end())
            return Integer(*p.second);
    }

    auto n(*this);

    if (n % 2)
        n += 2;
    else
        ++n;

    while (!n.isPrime())
        n += 2;

    return n;
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

bool Integer::isCoprime(Integer const& other) const noexcept
{
    return gcd(*this, other) == 1;
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
