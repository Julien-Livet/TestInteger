#include <QTest>

class TestInteger : public QObject
{
    Q_OBJECT

    private slots:
        void testAddition();
        void testAnd();
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
    QVERIFY(Integer<unsigned char>{0} + 0 == 0);
    QVERIFY(Integer<unsigned char>{0} + 1 == 1);
    QVERIFY(Integer<unsigned char>{1} + 2 == 3);
    QVERIFY(Integer<unsigned char>{1} + 3 == 4);
    QVERIFY(Integer<unsigned char>{-3} + 2 == -1);
    QVERIFY(Integer<unsigned char>{17} + 10 == 27);

    std::random_device rd;
    long long const a{rd()};
    long long const b{rd()};
    //std::cout << "a b a+b " << a << " " << b << " " << a + b << std::endl;
    QVERIFY((Integer<unsigned char>{a} + b).template cast<long long>() == a + b);
}

void TestInteger::testAnd()
{
    QVERIFY((Integer<unsigned char>{0} & 0) == 0);
    QVERIFY((Integer<unsigned char>{1} & 0) == 0);
    QVERIFY((Integer<unsigned char>{0} & 1) == 0);
    QVERIFY((Integer<unsigned char>{1} & 1) == 1);

    std::random_device rd;
    unsigned long long const a{rd()};
    unsigned long long const b{rd()};
    //std::cout << "a b a&b " << a << " " << b << " " << a & b << std::endl;
    QVERIFY((Integer<unsigned char>{a} & b).template cast<unsigned long long>() == (a & b));
}

void TestInteger::testChar()
{
    QVERIFY(!Integer<unsigned char>{'\0'});
    QVERIFY(Integer<unsigned char>{'\0'} == '\0');
    QVERIFY(Integer<unsigned char>{'1'} == '1');
    QVERIFY(Integer<unsigned char>{'a'}.cast<char>() == 'a');
}

void TestInteger::testDivision()
{
    QVERIFY(Integer<unsigned char>{2} / 3 == 0);
    QVERIFY(Integer<unsigned char>{3} / 2 == 1);
    QVERIFY(Integer<unsigned char>{4} / 2 == 2);
    QVERIFY(Integer<unsigned char>{-3} / 2 == -1);
    QVERIFY(Integer<unsigned char>{-4} / 2 == -2);
    QVERIFY(Integer<unsigned char>{5} / 2 == 2);
    QVERIFY(Integer<unsigned char>{6} / 2 == 3);
    QVERIFY(Integer<unsigned char>{6} / 3 == 2);

    {
        std::random_device rd;
        long long const a{rd()};
        long long const b{rd()};
        //std::cout << "a b a/b " << a << " " << b << " " << a / b << std::endl;
        QVERIFY((Integer<unsigned char>{a} / b).template cast<long long>() == a / b);
    }

    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        unsigned long const c{rd()};
        mpz_class const n{(mpz_class{a} * mpz_class{b}) / c};
        QVERIFY(((Integer<unsigned char>{a} * b) / c).toString() == n.get_str());
    }
}

void TestInteger::testEqualities()
{
    QVERIFY(-Integer<unsigned char>{1} == -1);
    QVERIFY(-Integer<unsigned long>{1} == -1);
    QVERIFY(Integer<unsigned char>{6} == 6);
}

void TestInteger::testInequalities()
{
    QVERIFY(Integer<unsigned char>{(short) 0} >= Integer<unsigned char>{(short) 0});
    QVERIFY(Integer<unsigned char>{(short) 0} <= Integer<unsigned char>{(short) 0});
    QVERIFY(Integer<unsigned char>{(short) 1} > Integer<unsigned char>{(short) 0});
    QVERIFY(Integer<unsigned char>{(short) -1} < Integer<unsigned char>{(short) 0});
    QVERIFY(Integer<unsigned char>{15} < 16);
    QVERIFY(Integer<unsigned char>{16} > 15);
    QVERIFY(Integer<unsigned char>{0} < 1);
    QVERIFY(Integer<unsigned char>{1} > 0);
    QVERIFY(Integer<unsigned char>{1} < 2);
    QVERIFY(Integer<unsigned char>{3} > 2);
    QVERIFY(Integer<unsigned char>{3} > 1);
    QVERIFY(Integer<unsigned char>{-1} < 0);
    QVERIFY(Integer<unsigned char>{-2} > -3);
}

void TestInteger::testModulo()
{
    QVERIFY(Integer<unsigned char>{2} % 2 == 0);
    QVERIFY(Integer<unsigned char>{4} % 2 == 0);
    QVERIFY(Integer<unsigned char>{5} % 3 == 2);
    QVERIFY(Integer<unsigned char>{-5} % -3 == -2);
    QVERIFY(Integer<unsigned char>{-5} % 3 == 1);
    QVERIFY(Integer<unsigned char>{5} % -3 == -1);
    QVERIFY(Integer<unsigned char>{5} % 1 == 0);

    {
        std::random_device rd;
        long long const a{rd()};
        long long const b{rd()};
        //std::cout << "a b a%b " << a << " " << b << " " << a % b << std::endl;
        QVERIFY((Integer<unsigned char>{a} % b).template cast<long long>() == a % b);
    }

    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        unsigned long const c{rd()};
        mpz_class const n{(mpz_class{a} * mpz_class{b}) % c};
        QVERIFY(((Integer<unsigned char>{a} * b) % c).toString() == n.get_str());
    }
}

void TestInteger::testMultiplication()
{
    QVERIFY(Integer<unsigned char>{0} * 1 == 0);
    QVERIFY(Integer<unsigned char>{0} * -1 == 0);
    QVERIFY(Integer<unsigned char>{1} * 0 == 0);
    QVERIFY(Integer<unsigned char>{-1} * 0 == 0);
    QVERIFY(Integer<unsigned char>{1} * 1 == 1);
    QVERIFY(Integer<unsigned char>{-1} * 1 == -1);
    QVERIFY(Integer<unsigned char>{-1} * -1 == 1);
    QVERIFY(Integer<unsigned char>{5} * 7 == 35);
    QVERIFY(Integer<unsigned char>{-5} * 7 == -35);
    QVERIFY(Integer<unsigned char>{5} * -7 == -35);
    QVERIFY(Integer<unsigned char>{-5} * -7 == 35);
    QVERIFY(-Integer<unsigned char>{5} * -7 == 35);
    QVERIFY(Integer<unsigned char>{827382986} * 2670051752 == 2209155391344291472ull);

    {
        std::random_device rd;
        long long const a{static_cast<short>(rd())};
        long long const b{static_cast<short>(rd())};
        //std::cout << "a b a*b " << a << " " << b << " " << a * b << std::endl;
        QVERIFY((Integer<unsigned char>{a} * b).template cast<long long>() == a * b);
    }

    {
        std::random_device rd;
        unsigned long const a{rd()};
        unsigned long const b{rd()};
        mpz_class const n{mpz_class{a} * mpz_class{b}};
        QVERIFY((Integer<unsigned char>{a} * b).toString() == n.get_str());
    }
}

void TestInteger::testOr()
{
    QVERIFY((Integer<unsigned char>{0} | 0) == 0);
    QVERIFY((Integer<unsigned char>{1} | 0) == 1);
    QVERIFY((Integer<unsigned char>{0} | 1) == 1);
    QVERIFY((Integer<unsigned char>{1} | 1) == 1);

    std::random_device rd;
    unsigned long long const a{rd()};
    unsigned long long const b{rd()};
    //std::cout << "a b a|b " << a << " " << b << " " << a | b << std::endl;
    QVERIFY((Integer<unsigned char>{a} | b).template cast<unsigned long long>() == (a | b));
}

void TestInteger::testPow()
{
    QVERIFY(pow(Integer<unsigned char>{0}, 0) == 1);
    QVERIFY(pow(Integer<unsigned char>{1}, 0) == 1);
    QVERIFY(pow(Integer<unsigned char>{2}, 0) == 1);
    QVERIFY(pow(Integer<unsigned char>{0}, 1) == 0);
    QVERIFY(pow(Integer<unsigned char>{0}, 2) == 0);
    QVERIFY(pow(Integer<unsigned char>{1}, 0) == 1);
    QVERIFY(pow(Integer<unsigned char>{1}, 1) == 1);
    QVERIFY(pow(Integer<unsigned char>{1}, 2) == 1);
    QVERIFY(pow(Integer<unsigned char>{2}, 2) == 4);
    QVERIFY(pow(Integer<unsigned char>{2}, 3) == 8);
    QVERIFY(pow(Integer<unsigned char>{3}, 2) == 9);
    QVERIFY(powm(Integer<unsigned char>{5}, 13, 23) == 21);
}

void TestInteger::testPrimes()
{
    QVERIFY(Integer<unsigned char>{(unsigned char)2}.isPrime());
    QVERIFY(!Integer<unsigned char>{(unsigned char)4}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)3}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)5}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)7}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)11}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)13}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)17}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)19}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)23}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)29}.isPrime());
    QVERIFY(Integer<unsigned char>{(unsigned char)31}.isPrime());
}

void TestInteger::testShift()
{
    QVERIFY((Integer<unsigned char>{1} >> 1) == 0);
    QVERIFY((Integer<unsigned char>{3} >> 1) == 1);
    QVERIFY((Integer<unsigned char>{1} << 1) == 2);
    QVERIFY((Integer<unsigned char>{1878050631} << 31) == 4033083020188581888ull);
    QVERIFY((Integer<unsigned char>{1800090071576084480ull} >> 21) == 858349834240ull);

    std::random_device rd;

    {
        unsigned long long const a{rd()};
        unsigned long long const b{rd() % 32};
        //std::cout << "a b a<<b " << a << " " << b << " " << (a << b) << std::endl;
        QVERIFY((Integer<unsigned char>{a} << b).template cast<long long>() == a << b);
    }

    {
        unsigned long long a{rd()};
        a <<= 32;
        unsigned long long const b{rd() % 32};
        //std::cout << "a b a>>b " << a << " " << b << " " << (a >> b) << std::endl;
        QVERIFY((Integer<unsigned char>{a} >> b).template cast<long long>() == a >> b);
    }
}

void TestInteger::testShort()
{
    QVERIFY(!Integer<unsigned char>{(short)0});
    QVERIFY(Integer<unsigned char>{(short)0} == (short)0);
    QVERIFY(Integer<unsigned char>{(short)1} == (short)1);
    QVERIFY(Integer<unsigned char>{(short)-1} == (short)-1);
    QVERIFY(Integer<unsigned char>{(short)0}.cast<short>() == (short)0);
    QVERIFY(Integer<unsigned char>{(short)1}.cast<short>() == (short)1);
    QVERIFY(Integer<unsigned char>{(short)-1}.cast<short>() == (short)-1);
}

void TestInteger::testString()
{
    QVERIFY((Integer<unsigned char>{"0b10", 2} == 2));
    QVERIFY((Integer<unsigned char>{"17", 10} == 17));
    QVERIFY((Integer<unsigned char>{"f", 16} == 15));
    QVERIFY((Integer<unsigned char>{"Z", 62} == 61));
}

void TestInteger::testSubstraction()
{
    QVERIFY(Integer<unsigned char>{-1} - 2 == -3);
    QVERIFY(Integer<unsigned char>{-1} - 3 == -4);
    QVERIFY(Integer<unsigned char>{-2} - 4 == -6);
    QVERIFY(Integer<unsigned char>{0} - 1 == -1);
    QVERIFY(Integer<unsigned char>{1} - 1 == 0);
    QVERIFY(Integer<unsigned char>{3} - 1 == 2);
    QVERIFY(Integer<unsigned char>{1} - 4 == -3);
    QVERIFY(Integer<unsigned char>{-5} - -3 == -2);
    QVERIFY(Integer<unsigned char>{3} - 2 == 1);
    QVERIFY(Integer<unsigned char>{17} - 10 == 7);

    std::random_device rd;
    long long const a{rd()};
    long long const b{rd()};
    //std::cout << "a b a-b " << a << " " << b << " " << a - b << std::endl;
    QVERIFY((Integer<unsigned char>{a} - b).template cast<long long>() == a - b);
}

void TestInteger::testToString()
{
    QVERIFY(Integer<unsigned char>{'\2'}.toString(2) == "0b00000010");
    QVERIFY(Integer<unsigned char>{(int)2}.toString(2) == "0b00000000000000000000000000000010");
    QVERIFY((Integer<unsigned char>{'\2'} >> 1).toString(2) == "0b00000001");
    QVERIFY((Integer<unsigned char>{'\2'} >> 2).toString(2) == "0b00000000");
    QVERIFY((Integer<unsigned char>{(unsigned char)255} << 1).toString(2) == "0b0000000111111110");
    QVERIFY(((Integer<unsigned char>{(unsigned char)255} << 1) >> 2).toString(2) == "0b0000000001111111");
    QVERIFY(Integer<unsigned char>{10}.toString() == "10");
    QVERIFY(Integer<unsigned char>{17}.toString() == "17");
    QVERIFY(Integer<unsigned char>{15}.toString(16) == "0xf");
    QVERIFY(Integer<unsigned char>{10}.toString(62) == "a");
    QVERIFY(Integer<unsigned char>{35}.toString(62) == "z");
    QVERIFY(Integer<unsigned char>{36}.toString(62) == "A");
    QVERIFY(Integer<unsigned char>{61}.toString(62) == "Z");
    QVERIFY(Integerll{"94882770404074832568499"}.toString() == mpz_class{"94882770404074832568499"}.get_str());
}

void TestInteger::testUnsignedLongLong()
{
    QVERIFY(!Integer<unsigned long long>{0ull});
    QVERIFY(Integer<unsigned long long>{1ull} == 1);
    QVERIFY(Integer<unsigned long long>{1ull}.cast<unsigned long long>() == 1ull);
}

void TestInteger::testUnsignedShort()
{
    QVERIFY(!Integer<unsigned char>{(unsigned short)0});
    QVERIFY(Integer<unsigned char>{(unsigned short)0} == (unsigned short)0);
    QVERIFY(Integer<unsigned char>{(unsigned short)1} == (unsigned short)1);
    QVERIFY(Integer<unsigned char>{(unsigned short)-1} == (unsigned short)-1);
    QVERIFY(Integer<unsigned char>{(unsigned short)0}.cast<unsigned short>() == (unsigned short)0);
    QVERIFY(Integer<unsigned char>{(unsigned short)1}.cast<unsigned short>() == (unsigned short)1);
    QVERIFY(Integer<unsigned char>{(unsigned short)-1}.cast<unsigned short>() == (unsigned short)-1);
}

QTEST_MAIN(TestInteger)
#include "testinteger.moc"
