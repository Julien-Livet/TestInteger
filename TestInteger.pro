QT += testlib

CONFIG += c++17 object_parallel_to_source

CONFIG(debug, debug|release) {
} else {
    DEFINES += NDEBUG
}

DEFINES += USING_GMP

INCLUDEPATH += .. \
            D:/Programmes/gmp-6.3.0

SOURCES = testinteger.cpp

HEADERS = ../Integer.h

LIBS += D:/Programmes/gmp-6.3.0/.libs/libgmp.dll.a
