#include <QTest>

class TestInteger : public QObject
{
    Q_OBJECT

    private slots:
        void testAddition();
};

#include <iostream>

#include "Integer.h"

void TestInteger::testAddition()
{
    QVERIFY(Integer(0) + 0 == 0);
}

QTEST_MAIN(TestInteger)
#include "testinteger.moc"
