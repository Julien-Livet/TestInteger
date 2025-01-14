#include <QTest>

class TestInteger : public QObject
{
    Q_OBJECT

    private slots:
        void testAddition();
        void testAnd();
        void testBits();
        void testChar();
        void testDivision();
        void testEqualities();
        void testGcd();
        void testInequalities();
        void testModulo();
        void testMultiplication();
        void testOr();
        void testPow();
        void testPrimes();
        void testShift();
        void testShort();
        void testString();
        void testSubstraction();
        void testToString();
        void testUnsignedLongLong();
        void testUnsignedShort();
};

#include <iostream>

#include "Integer.h"

void TestInteger::testAddition()
{
    QVERIFY(Integerc(0) + 0 == 0);
    QVERIFY(Integerc(0) + 1 == 1);
    QVERIFY(Integerc(1) + 2 == 3);
    QVERIFY(Integerc(1) + 3 == 4);
    QVERIFY(Integerc(-3) + 2 == -1);
    QVERIFY(Integerc(17) + 10 == 27);

#ifdef WITH_GMP
    {
        Integerll n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom<std::random_device>();
        n2.setRandom<std::random_device>();

        auto const n3(n1 + n2);
        auto const n1_{n1.cast<mpz_class>()};
        auto const n2_{n2.cast<mpz_class>()};
        mpz_class const n3_{n1_ + n2_};

        QVERIFY(n3 == n3_.get_str());
    }
#endif

    std::random_device rd;
    long long const a{rd()};
    long long const b{rd()};
    //std::cout << "a b a+b " << a << " " << b << " " << a + b << std::endl;
    QVERIFY((Integerc(a) + b).cast<long long>() == a + b);
}

void TestInteger::testAnd()
{
    QVERIFY((Integerc(0) & 0) == 0);
    QVERIFY((Integerc(1) & 0) == 0);
    QVERIFY((Integerc(0) & 1) == 0);
    QVERIFY((Integerc(1) & 1) == 1);

#ifdef WITH_GMP
    {
        Integerll n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom<std::random_device>();
        n1.setPositive();
        n2.setRandom<std::random_device>();
        n2.setPositive();

        auto const n3(n1 & n2);
        auto const n1_{n1.cast<mpz_class>()};
        auto const n2_{n2.cast<mpz_class>()};
        mpz_class const n3_{n1_ & n2_};

        QVERIFY(n3 == n3_.get_str());
    }
#endif

    std::random_device rd;
    unsigned long long const a{rd()};
    unsigned long long const b{rd()};
    //std::cout << "a b a&b " << a << " " << b << " " << a & b << std::endl;
    QVERIFY((Integerc(a) & b).cast<unsigned long long>() == (a & b));
}

void TestInteger::testBits()
{
    {
        std::string const s{"0b10001000101011110010"};
        Integerc const n(s, 2);

        size_t k{0};

        while (s[s.size() - 1 - k] != 'b')
        {
            QVERIFY(static_cast<char>('0' + n.bit(k)) == s[s.size() - 1 - k]);
            ++k;
        }

        QVERIFY(n.count() == 9);
        QVERIFY(n.number() == 20);
    }

    {
        std::string const s{"0b10001000101011110010"};
        Integerc n('\0');

        size_t k{0};

        while (s[s.size() - 1 - k] != 'b')
        {
            n.setBit(k, s[s.size() - 1 - k] == '1');
            ++k;
        }

        QVERIFY((n == Integerc{s, 2}));
    }
}

void TestInteger::testChar()
{
    QVERIFY(!Integerc('\0'));
    QVERIFY(Integerc('\0') == '\0');
    QVERIFY(Integerc('1') == '1');
    QVERIFY(Integerc('a').cast<char>() == 'a');
}

void TestInteger::testDivision()
{
    QVERIFY(Integerc(12) / 4 == 3);
    QVERIFY(Integerc(2) / 3 == 0);
    QVERIFY(Integerc(3) / 2 == 1);
    QVERIFY(Integerc(4) / 2 == 2);
    QVERIFY(Integerc(-3) / 2 == -1);
    QVERIFY(Integerc(-4) / 2 == -2);
    QVERIFY(Integerc(5) / 2 == 2);
    QVERIFY(Integerc(6) / 2 == 3);
    QVERIFY(Integerc(6) / 3 == 2);
    QVERIFY((Integerc{122, 17, 200,43} / Integerc{23, 117} == Integerc{5, 52, 54}));

    {
        Integer32 const a{0, 561453, 13205, 1564};
        Integer32 const b{1698, 721};

        //auto const qr{computeQrBinary(a, b)};
        auto const qr{computeQrBurnikelZiegler(a, b)};

        QVERIFY((qr.first == Integer32{330, 2815252282}));
        QVERIFY((qr.second == Integer32{414, 1722637250}));
    }

    {
        Integerll const a("0b10101011110100110001100110101111100010110000011001001100100101110100111011100110100100101110001010001100110111011001101010011011000010111111001010110101010");
        Integerll const b("0b100010001011011011011011001110101100100110101011101100101000011011111110010011010101101010101011111110010101010");

        auto const qr{computeQrBurnikelZiegler(a, b)};

        QVERIFY((qr.first == Integerll{22110129672729ull}));
        QVERIFY((qr.second == Integerll{19184386057769ull, 9219185875676460304ull}));
    }

    {
        Integerc const a{122, 17, 200,43};
        Integerc const b{23, 117};

        auto const qr{computeQrBurnikelZiegler(a, b)};

        QVERIFY((qr.first == Integerc{5, 52, 54}));
        QVERIFY((qr.second == Integerc{17, 125}));
    }

    {
        auto const qr1{computeQr(Integerc(10), 3)};
        QVERIFY(qr1.first == 3);
        QVERIFY(qr1.second == 1);

        auto const qr2{computeQr(Integerc(10), -3)};
        QVERIFY(qr2.first == -4);
        QVERIFY(qr2.second == -2);

        auto const qr3{computeQr(Integerc(-10), -3)};
        QVERIFY(qr3.first == 4);
        QVERIFY(qr3.second == -2);

        auto const qr4{computeQr(Integerc(-10), 3)};
        QVERIFY(qr4.first == -4);
        QVERIFY(qr4.second == 2);
    }

#ifdef WITH_GMP
    {
        Integerll n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom<std::random_device>();
        n1.setPositive();

        do
        {
            n2.setRandom<std::random_device>();
            n2.setPositive();
        } while (!n2);

        auto const q(n1 / n2);
        auto const r(n1 % n2);
        auto const n1_{n1.cast<mpz_class>()};
        auto const n2_{n2.cast<mpz_class>()};
        mpz_class q_, r_;

        mpz_divmod(q_.get_mpz_t(), r_.get_mpz_t(), n1_.get_mpz_t(), n2_.get_mpz_t());

        QVERIFY(q == q_.get_str());
        QVERIFY(r == r_.get_str());
        QVERIFY(n1 == q * n2 + r);
        QVERIFY(n1_ == q_ * n2_ + r_);
    }
#endif

    {
        std::random_device rd;
        long long const a{rd()};
        long long const b{rd()};
        //std::cout << "a b a/b " << a << " " << b << " " << a / b << std::endl;
        QVERIFY((Integerc(a) / b).cast<long long>() == a / b);
    }

#ifdef WITH_GMP
    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        unsigned long const c{rd()};
        mpz_class const n{(mpz_class{a} * mpz_class{b}) / c};
        QVERIFY((Integerc(a) * b) / c == n.get_str());
    }
#endif
}

void TestInteger::testEqualities()
{
    QVERIFY(-Integerc(1) == -1);
    QVERIFY(-Integerl(1) == -1);
    QVERIFY(Integerc(6) == 6);

    Integerll const a(234);
    Integerll const b(167);
    Integerll const n(293);

    auto reduction{[] (Integerll const& t, Integerll const& R, Integerll const& n, Integerll const& n_) -> Integerll
        {
            auto const m(((t % R) * n_) % R);
            auto const x((t + m * n) / R);

            if (x < n)
                return x;
            else
                return x - n;
        }
    };

    auto redmulmod{[&reduction] (Integerll const& a, Integerll const& b, Integerll const& n,
                                 Integerll const& R, Integerll const& n_, Integerll const& R2modn) -> Integerll
        {
            auto const reda(reduction(a * R2modn, R, n, n_));
            auto const redb(reduction(b * R2modn, R, n, n_));
            auto const redc(reduction(reda * redb, R, n, n_));

            return reduction(redc, R, n, n_);
        }
    };

    {
        Integerll const R(1000);
        Integerll const R_(247), n_(843);

        auto const R2modn((R * R) % n);

        QVERIFY(R2modn == 284);

        QVERIFY(((a * b) % n == redmulmod(a, b, n, R, n_, R2modn)));
    }

    {
        Integerll R(2);

        while (R <= n)
            R <<= 1;

        if (!(n & 1))
            ++R;

        while (!n.isCoprime(R))
        {
            if (!(n & 1))
                --R;

            R <<= 1;

            if (!(n & 1))
                ++R;
        }

        Integerll R_, n_;

        auto const d(gcdExtended(R, -n, R_, n_));

        QVERIFY(R * R_ - n * n_ == d);

        if (d == -1)
        {
            R_ = -R;
            n_ = -n_;
        }

        auto const R2modn((R * R) % n);

        QVERIFY(((a * b) % n == redmulmod(a, b, n, R, n_, R2modn)));
    }

    {
        Integerll x(4), y(16), m(23);
        Integerll R(2);

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

        Integerll R_, m_;

        auto const d(gcdExtended(R, -m, R_, m_));

        QVERIFY(R * R_ - m * m_ == d);

        if (d == -1)
        {
            R_ = -R;
            m_ = -m_;
        }

        auto const R2modm((R * R) % m);

        auto res(redmulmod(x, y, m, R, m_, R2modm));
        while (res < 0)
            res += m;

        QVERIFY((x * y) % m == res);
    }
}

void TestInteger::testGcd()
{
    QVERIFY(gcd(8_z, 12) == 4);

    auto const a{120_z};
    auto const b{23_z};
    auto u{0_z}, v{0_z};

    auto const d(gcdExtended(a, b, u, v));

    QVERIFY(d == 1);
    QVERIFY(a * u + b * v == d);

    QVERIFY(lcm(4_z, 6) == 12);
    QVERIFY(lcm(21_z, 6) == 42);
}

void TestInteger::testInequalities()
{
    QVERIFY(Integerc((short)0) >= Integerc{(short)0});
    QVERIFY(Integerc((short)0) <= Integerc{(short)0});
    QVERIFY(Integerc((short)1) > Integerc{(short)0});
    QVERIFY(Integerc((short)-1) < Integerc{(short)0});
    QVERIFY(Integerc(15) < 16);
    QVERIFY(Integerc(16) > 15);
    QVERIFY(Integerc(0) < 1);
    QVERIFY(Integerc(1) > 0);
    QVERIFY(Integerc(1) < 2);
    QVERIFY(Integerc(3) > 2);
    QVERIFY(Integerc(3) > 1);
    QVERIFY(Integerc(-1) < 0);
    QVERIFY(Integerc(-2) > -3);
}

void TestInteger::testModulo()
{
    QVERIFY(Integerc(2) % 2 == 0);
    QVERIFY(Integerc(4) % 2 == 0);
    QVERIFY(Integerc(5) % 3 == 2);
    QVERIFY(Integerc(-5) % -3 == -2);
    QVERIFY(Integerc(-5) % 3 == 1);
    QVERIFY(Integerc(5) % -3 == -1);
    QVERIFY(Integerc(5) % 1 == 0);

#ifdef WITH_GMP
    {
        Integerll n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom<std::random_device>();
        n1.setPositive();

        do
        {
            n2.setRandom<std::random_device>();
            n2.setPositive();
        } while (!n2);

        auto const n3(n1 % n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class const n3_{n1_ % n2_};

        QVERIFY(n3 == n3_.get_str());
    }
#endif

    {
        std::random_device rd;
        long long const a{rd()};
        long long const b{rd()};
        //std::cout << "a b a%b " << a << " " << b << " " << a % b << std::endl;
        QVERIFY((Integerc(a) % b).cast<long long>() == a % b);
    }

#ifdef WITH_GMP
    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        unsigned long const c{rd()};
        mpz_class const n{(mpz_class{a} * mpz_class{b}) % c};
        QVERIFY((Integerc(a) * b) % c == n.get_str());
    }
#endif
}

void TestInteger::testMultiplication()
{
    QVERIFY(Integerc(0) * 1 == 0);
    QVERIFY(Integerc(0) * -1 == 0);
    QVERIFY(Integerc(1) * 0 == 0);
    QVERIFY(Integerc(-1) * 0 == 0);
    QVERIFY(Integerc(1) * 1 == 1);
    QVERIFY(Integerc(-1) * 1 == -1);
    QVERIFY(Integerc(-1) * -1 == 1);
    QVERIFY(Integerc(5) * 7 == 35);
    QVERIFY(Integerc(-5) * 7 == -35);
    QVERIFY(Integerc(5) * -7 == -35);
    QVERIFY(Integerc(-5) * -7 == 35);
    QVERIFY(-Integerc(5) * -7 == 35);
    QVERIFY(Integerc(827382986) * 2670051752 == 2209155391344291472ull);
    QVERIFY(Integerll(17539966127645434034ull) * 452240028370ull == Integerll("7932274779175210128837123544580"));

#ifdef WITH_GMP
    {
        Integerll n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom<std::random_device>();
        n2.setRandom<std::random_device>();

        auto const n3(n1 * n2);
        auto const n1_{n1.cast<mpz_class>()};
        auto const n2_{n2.cast<mpz_class>()};
        mpz_class const n3_{n1_ * n2_};

        QVERIFY(n3 == n3_.get_str());
    }
#endif

    {
        std::random_device rd;
        long long const a{static_cast<short>(rd())};
        long long const b{static_cast<short>(rd())};
        //std::cout << "a b a*b " << a << " " << b << " " << a * b << std::endl;
        QVERIFY((Integerc(a) * b).cast<long long>() == a * b);
    }

#ifdef WITH_GMP
    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        mpz_class const n{mpz_class{a} * mpz_class{b}};
        QVERIFY(Integerc(a) * b == n.get_str());
    }

    {
        Integerll const a{0ull, 0ull, 35521434ull, 14919252733983618111ull};
        Integerll const b{11163967620057660708ull, 5245663108828415198ull, 4063525853430116991ull, 8147888677444386816ull};

        QVERIFY(a * b == "45918528047859382727956397566465766535802613706550573885714125513512139783633707483925386458498848522240");
        QVERIFY(a * b == mpz_class{a.cast<mpz_class>() * b.cast<mpz_class>()});
    }

    {
        Integer32 const a{330, 2815252282};
        Integer32 const b{1698, 721};

        QVERIFY(a * b == mpz_class{a.cast<mpz_class>() * b.cast<mpz_class>()});
    }
#endif
}

void TestInteger::testOr()
{
    QVERIFY((Integerc(0) | 0) == 0);
    QVERIFY((Integerc(1) | 0) == 1);
    QVERIFY((Integerc(0) | 1) == 1);
    QVERIFY((Integerc(1) | 1) == 1);

#ifdef WITH_GMP
    {
        Integerll n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom<std::random_device>();
        n1.setPositive();
        n2.setRandom<std::random_device>();
        n2.setPositive();

        auto const n3(n1 | n2);
        auto const n1_{n1.cast<mpz_class>()};
        auto const n2_{n2.cast<mpz_class>()};
        mpz_class const n3_{n1_ | n2_};

        QVERIFY(n3 == n3_.get_str());
    }
#endif

    std::random_device rd;
    unsigned long long const a{rd()};
    unsigned long long const b{rd()};
    //std::cout << "a b a|b " << a << " " << b << " " << a | b << std::endl;
    QVERIFY((Integerc(a) | b).cast<unsigned long long>() == (a | b));
}

void TestInteger::testPow()
{
    QVERIFY(pow(Integerc(0), 0) == 1);
    QVERIFY(pow(Integerc(1), 0) == 1);
    QVERIFY(pow(Integerc(2), 0) == 1);
    QVERIFY(pow(Integerc(0), 1) == 0);
    QVERIFY(pow(Integerc(0), 2) == 0);
    QVERIFY(pow(Integerc(1), 0) == 1);
    QVERIFY(pow(Integerc(1), 1) == 1);
    QVERIFY(pow(Integerc(1), 2) == 1);
    QVERIFY(pow(Integerc(2), 2) == 4);
    QVERIFY(pow(Integerc(2), 3) == 8);
    QVERIFY(pow(Integerc(3), 2) == 9);
    QVERIFY(powm(Integerc(5), 13, 23) == 21);
}

void TestInteger::testPrimes()
{
#ifdef WITH_GMP
    {
        mpz_class n;

        mpz_primorial_ui(n.get_mpz_t(), 0);
        QVERIFY(primorial(0_z) == n);
        mpz_primorial_ui(n.get_mpz_t(), 1);
        QVERIFY(primorial(1_z) == n);
        mpz_primorial_ui(n.get_mpz_t(), 2);
        QVERIFY(primorial(2_z) == n);
        mpz_primorial_ui(n.get_mpz_t(), 7);
        QVERIFY(primorial(7_z) == n);
    }

    {
        mpz_class n1{2}, n2{9};

        QVERIFY(mpz_jacobi(n1.get_mpz_t(), n2.get_mpz_t()) == jacobi(Integerll(n1), n2));
    }
#endif

    QVERIFY((2_z).isPrime() == 2);
    QVERIFY(!(4_z).isPrime());
    QVERIFY((3_z).isPrime() == 2);
    QVERIFY((3_z).nextPrime() == 5);
    QVERIFY((4_z).nextPrime() == 5);
    QVERIFY((5_z).nextPrime() == 7);
    QVERIFY((10_z).nextPrime() == 11);
    QVERIFY((13_z).previousPrime() == 11);
    QVERIFY((12_z).previousPrime() == 11);
    QVERIFY((11_z).previousPrime() == 7);
    QVERIFY((5_z).isPrime() == 2);
    QVERIFY((7_z).isPrime() == 2);
    QVERIFY((11_z).isPrime() == 2);
    QVERIFY((13_z).isPrime() == 2);
    QVERIFY((17_z).isPrime() == 2);
    QVERIFY((19_z).isPrime() == 2);
    QVERIFY((23_z).isPrime() == 2);
    QVERIFY((29_z).isPrime() == 2);
    QVERIFY((31_z).isPrime() == 2);
    QVERIFY((40'562_z).nextPrime().isPrime() == 2);
    QVERIFY((15'465'319_z).previousPrime().isPrime() == 2);
    QVERIFY((15'465'319_z).nextPrime().isPrime() == 2);
    QVERIFY((1'299'709_z).isPrime());//100'000th prime
    QVERIFY(!(15'482'009_z).isPrime());
    QVERIFY((13'359'555'403_z).isPrime());//600'000'000th prime
    QVERIFY((4113101149215104800030529537915953170486139623539759933135949994882770404074832568499_z).isPrime());
}

void TestInteger::testShift()
{
    QVERIFY((Integerc(1) >> 1) == 0);
    QVERIFY((Integerc(3) >> 1) == 1);
    QVERIFY((Integerc(1) << 1) == 2);
    QVERIFY((Integerc(1878050631) << 31) == 4033083020188581888ull);
    QVERIFY((Integerc(1800090071576084480ull) >> 21) == 858349834240ull);

    std::random_device rd;

    {
        unsigned long long const a{rd()};
        unsigned long long const b{rd() % 32};
        //std::cout << "a b a<<b " << a << " " << b << " " << (a << b) << std::endl;
        QVERIFY((Integerc(a) << b).cast<unsigned long long>() == a << b);
    }

    {
        unsigned long long a{rd()};
        a <<= 32;
        unsigned long long const b{rd() % 32};
        //std::cout << "a b a>>b " << a << " " << b << " " << (a >> b) << std::endl;
        QVERIFY((Integerc(a) >> b).cast<unsigned long long>() == a >> b);
    }
}

void TestInteger::testShort()
{
    QVERIFY(!Integerc((short)0));
    QVERIFY(Integerc((short)0) == (short)0);
    QVERIFY(Integerc((short)1) == (short)1);
    QVERIFY(Integerc((short)-1) == (short)-1);
    QVERIFY(Integerc((short)0).cast<short>() == (short)0);
    QVERIFY(Integerc((short)1).cast<short>() == (short)1);
    QVERIFY(Integerc((short)-1).cast<short>() == (short)-1);
}

void TestInteger::testString()
{
    QVERIFY((Integerc("0b10", 2) == 2));
    QVERIFY((Integerc("17", 10) == 17));
    QVERIFY((Integerc("f", 16) == 15));
    QVERIFY((Integerc("Z", 62) == 61));
}

void TestInteger::testSubstraction()
{
    QVERIFY(Integerc(-1) - 2 == -3);
    QVERIFY(Integerc(-1) - 3 == -4);
    QVERIFY(Integerc(-2) - 4 == -6);
    QVERIFY(Integerc(0) - 1 == -1);
    QVERIFY(Integerc(1) - 1 == 0);
    QVERIFY(Integerc(3) - 1 == 2);
    QVERIFY(Integerc(1) - 4 == -3);
    QVERIFY(Integerc(-5) - -3 == -2);
    QVERIFY(Integerc(3) - 2 == 1);
    QVERIFY(Integerc(17) - 10 == 7);

#ifdef WITH_GMP
    {
        Integerll n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom<std::random_device>();
        n2.setRandom<std::random_device>();

        auto const n3(n1 - n2);
        auto const n1_{n1.cast<mpz_class>()};
        auto const n2_{n2.cast<mpz_class>()};
        mpz_class const n3_{n1_ - n2_};

        QVERIFY(n3 == n3_.get_str());
    }
#endif

    std::random_device rd;
    long long const a{rd()};
    long long const b{rd()};
    //std::cout << "a b a-b " << a << " " << b << " " << a - b << std::endl;
    QVERIFY((Integerc(a) - b).cast<long long>() == a - b);
}

void TestInteger::testToString()
{
    QVERIFY(Integerc('\2').toString(2) == "0b00000010");
    QVERIFY(Integerc((int)2).toString(2) == "0b00000010");
    QVERIFY((Integerc('\2') >> 1).toString(2) == "0b00000001");
    QVERIFY((Integerc('\2') >> 2).toString(2) == "0b00000000");
    QVERIFY((Integerc((unsigned char)255) << 1).toString(2) == "0b0000000111111110");
    QVERIFY(((Integerc((unsigned char)255) << 1) >> 2).toString(2) == "0b01111111");
    QVERIFY(Integerc(10).toString() == "10");
    QVERIFY(Integerc(17).toString() == "17");
    QVERIFY(Integerc(15).toString(16) == "0xf");
    QVERIFY(Integerc(10).toString(62) == "a");
    QVERIFY(Integerc(35).toString(62) == "z");
    QVERIFY(Integerc(36).toString(62) == "A");
    QVERIFY(Integerc(61).toString(62) == "Z");
    QVERIFY(Integerll("94882770404074832568499").toString() == "94882770404074832568499");
}

void TestInteger::testUnsignedLongLong()
{
    QVERIFY(!Integerll(0ull));
    QVERIFY(Integerll(1ull) == 1);
    QVERIFY(Integerll(1ull).cast<unsigned long long>() == 1ull);

    {
        Integerll n(~0ull);

        QVERIFY(n.bits().size() == 1 && n.bits().back() == ~0ull);

        n += 1;

        QVERIFY(n.bits().size() == 2 && n.bits()[0] == 1 && n.bits()[1] == 0);

        n -= 1;

        QVERIFY(n.bits().size() == 1 && n.bits().back() == ~0ull);
    }

    {
        Integerll const p(pow(Integerll(10), 20));

        QVERIFY(p.bits().size() == 2 && p.bits()[0] == 5 && p.bits()[1] == 7766279631452241920ull);
    }

    {
        Integerll p(pow(Integerll(10), 21));

        p *= 10;

        QVERIFY(p.bits().size() == 2 && p.bits()[0] == 542 && p.bits()[1] == 1864712049423024128ull);
    }

    {
        Integerll const dividend{16064214685684893896ull, 1233950369490649711ull, 6706106583984356886ull,
                                 13049457702161808613ull, 14797807619686341699ull};
        Integerll const divisor{35521434ull, 14919252733983618111ull, 1302913595559511957ull,
                                9115028167675518012ull, 17539966127645434034ull};

        QVERIFY(divisor * 452240028370ull == "1860108980409718597998329855093233652496243230041402151984391328752685442362998644380535488288260");

        auto const qr{computeQr(dividend, divisor)};

        QVERIFY(qr.first == 452240028370ull);
        QVERIFY(qr.second == "2511973644799747096503184818228555377241860558862408039904124816691570905857945011775");
    }
}

void TestInteger::testUnsignedShort()
{
    QVERIFY(!Integers((unsigned short)0));
    QVERIFY(Integers((unsigned short)0) == (unsigned short)0);
    QVERIFY(Integers((unsigned short)1) == (unsigned short)1);
    QVERIFY(Integers((unsigned short)-1) == (unsigned short)-1);
    QVERIFY(Integers((unsigned short)0).cast<unsigned short>() == (unsigned short)0);
    QVERIFY(Integers((unsigned short)1).cast<unsigned short>() == (unsigned short)1);
    QVERIFY(Integers((unsigned short)-1).cast<unsigned short>() == (unsigned short)-1);
}

QTEST_MAIN(TestInteger)
#include "testinteger.moc"
