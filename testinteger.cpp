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

#include <gmpxx.h>

#include "Integer.h"

void TestInteger::testAddition()
{
    QVERIFY(Integerc(0) + 0 == 0);
    QVERIFY(Integerc(0) + 1 == 1);
    QVERIFY(Integerc(1) + 2 == 3);
    QVERIFY(Integerc(1) + 3 == 4);
    QVERIFY(Integerc(-3) + 2 == -1);
    QVERIFY(Integerc(17) + 10 == 27);

    {
        std::random_device rd;
        Integerc n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom(rd);
        n1.setPositive();

        auto const n3(n1 + n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class const n3_{n1_ + n2_};

        QVERIFY(n3 == n3_.get_str());
    }

    std::random_device rd;
    long long const a{rd()};
    long long const b{rd()};
    //std::cout << "a b a+b " << a << " " << b << " " << a + b << std::endl;
    QVERIFY((Integerc(a) + b).template cast<long long>() == a + b);
}

void TestInteger::testAnd()
{
    QVERIFY((Integerc(0) & 0) == 0);
    QVERIFY((Integerc(1) & 0) == 0);
    QVERIFY((Integerc(0) & 1) == 0);
    QVERIFY((Integerc(1) & 1) == 1);

    {
        std::random_device rd;
        Integerc n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom(rd);
        n1.setPositive();

        auto const n3(n1 & n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class const n3_{n1_ & n2_};

        QVERIFY(n3 == n3_.get_str());
    }

    std::random_device rd;
    unsigned long long const a{rd()};
    unsigned long long const b{rd()};
    //std::cout << "a b a&b " << a << " " << b << " " << a & b << std::endl;
    QVERIFY((Integerc(a) & b).template cast<unsigned long long>() == (a & b));
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

    {
        std::random_device rd;
        Integerc n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom(rd);
        n1.setPositive();

        do
        {
            n2.setRandom(rd);
            n2.setPositive();
        } while (!n2);

        Integerc const q(n1 / n2);
        Integerc const r(n1 % n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class q_, r_;

        mpz_divmod(q_.get_mpz_t(), r_.get_mpz_t(), n1_.get_mpz_t(), n2_.get_mpz_t());

        QVERIFY(q == q_.get_str());
        QVERIFY(r == r_.get_str());
        QVERIFY(n1 == q * n2 + r);
        QVERIFY(n1_ == q_ * n2_ + r_);
    }

    {
        std::random_device rd;
        long long const a{rd()};
        long long const b{rd()};
        //std::cout << "a b a/b " << a << " " << b << " " << a / b << std::endl;
        QVERIFY((Integerc(a) / b).template cast<long long>() == a / b);
    }

    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        unsigned long const c{rd()};
        mpz_class const n{(mpz_class{a} * mpz_class{b}) / c};
        QVERIFY((Integerc(a) * b) / c == n.get_str());
    }
}

void TestInteger::testEqualities()
{
    QVERIFY(-Integerc(1) == -1);
    QVERIFY(-Integerl(1) == -1);
    QVERIFY(Integerc(6) == 6);
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

    {
        std::random_device rd;
        Integerc n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom(rd);
        n1.setPositive();

        do
        {
            n2.setRandom(rd);
            n2.setPositive();
        } while (!n2);

        auto const n3(n1 % n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class const n3_{n1_ % n2_};

        QVERIFY(n3 == n3_.get_str());
    }

    {
        std::random_device rd;
        long long const a{rd()};
        long long const b{rd()};
        //std::cout << "a b a%b " << a << " " << b << " " << a % b << std::endl;
        QVERIFY((Integerc(a) % b).template cast<long long>() == a % b);
    }

    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        unsigned long const c{rd()};
        mpz_class const n{(mpz_class{a} * mpz_class{b}) % c};
        QVERIFY((Integerc(a) * b) % c == n.get_str());
    }
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

    {
        std::random_device rd;
        Integerc n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom(rd);
        n1.setPositive();

        auto const n3(n1 * n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class const n3_{n1_ * n2_};

        QVERIFY(n3 == n3_.get_str());
    }

    {
        std::random_device rd;
        long long const a{static_cast<short>(rd())};
        long long const b{static_cast<short>(rd())};
        //std::cout << "a b a*b " << a << " " << b << " " << a * b << std::endl;
        QVERIFY((Integerc(a) * b).template cast<long long>() == a * b);
    }

    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        mpz_class const n{mpz_class{a} * mpz_class{b}};
        QVERIFY(Integerc(a) * b == n.get_str());
    }
}

void TestInteger::testOr()
{
    QVERIFY((Integerc(0) | 0) == 0);
    QVERIFY((Integerc(1) | 0) == 1);
    QVERIFY((Integerc(0) | 1) == 1);
    QVERIFY((Integerc(1) | 1) == 1);

    {
        std::random_device rd;
        Integerc n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom(rd);
        n1.setPositive();

        auto const n3(n1 | n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class const n3_{n1_ | n2_};

        QVERIFY(n3 == n3_.get_str());
    }

    std::random_device rd;
    unsigned long long const a{rd()};
    unsigned long long const b{rd()};
    //std::cout << "a b a|b " << a << " " << b << " " << a | b << std::endl;
    QVERIFY((Integerc(a) | b).template cast<unsigned long long>() == (a | b));
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
    QVERIFY(Integerc((unsigned char)2).isPrime());
    QVERIFY(!Integerc((unsigned char)4).isPrime());
    QVERIFY(Integerc((unsigned char)3).isPrime());
    QVERIFY(Integerc((unsigned char)5).isPrime());
    QVERIFY(Integerc((unsigned char)7).isPrime());
    QVERIFY(Integerc((unsigned char)11).isPrime());
    QVERIFY(Integerc((unsigned char)13).isPrime());
    QVERIFY(Integerc((unsigned char)17).isPrime());
    QVERIFY(Integerc((unsigned char)19).isPrime());
    QVERIFY(Integerc((unsigned char)23).isPrime());
    QVERIFY(Integerc((unsigned char)29).isPrime());
    QVERIFY(Integerc((unsigned char)31).isPrime());
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
        QVERIFY((Integerc(a) << b).template cast<unsigned long long>() == a << b);
    }

    {
        unsigned long long a{rd()};
        a <<= 32;
        unsigned long long const b{rd() % 32};
        //std::cout << "a b a>>b " << a << " " << b << " " << (a >> b) << std::endl;
        QVERIFY((Integerc(a) >> b).template cast<unsigned long long>() == a >> b);
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

    {
        std::random_device rd;
        Integerc n1, n2;
        n1.setPrecision(4);
        n2.setPrecision(4);
        n1.setRandom(rd);
        n1.setPositive();

        auto const n3(n1 - n2);
        mpz_class const n1_{n1.toString(2).substr(2), 2};
        mpz_class const n2_{n2.toString(2).substr(2), 2};
        mpz_class const n3_{n1_ - n2_};

        QVERIFY(n3 == n3_.get_str());
    }

    std::random_device rd;
    long long const a{rd()};
    long long const b{rd()};
    //std::cout << "a b a-b " << a << " " << b << " " << a - b << std::endl;
    QVERIFY((Integerc(a) - b).template cast<long long>() == a - b);
}

void TestInteger::testToString()
{
    QVERIFY(Integerc('\2').toString(2) == "0b00000010");
    QVERIFY(Integerc((int)2).toString(2) == "0b00000000000000000000000000000010");
    QVERIFY((Integerc('\2') >> 1).toString(2) == "0b00000001");
    QVERIFY((Integerc('\2') >> 2).toString(2) == "0b00000000");
    QVERIFY((Integerc((unsigned char)255) << 1).toString(2) == "0b0000000111111110");
    QVERIFY(((Integerc((unsigned char)255) << 1) >> 2).toString(2) == "0b0000000001111111");
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

        QVERIFY(n.bits().size() == 2 && n.bits()[0] == 0 && n.bits()[1] == ~0ull);
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

        auto const qr{computeQr(dividend, divisor)};

        std::cout << "q " << qr.first << std::endl;
        std::cout << "r " << qr.second << std::endl;
        //QVERIFY(qr.first == 452240028370ull);
        //QVERIFY(qr.second == "2511973644799747096503184818228555377241860558862408039904124816691570905857945011775");
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
