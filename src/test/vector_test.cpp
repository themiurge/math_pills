#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../core/core.hpp"

TEST_CASE( "Empty vector creation", "[vector]" ) {
    mp::vector v;
    REQUIRE( v.getRows() == 0 );
    REQUIRE( v.getCols() == 1 );
}

TEST_CASE( "Vector append", "[vector]" ) {
    mp::vector v;
    REQUIRE( v.getRows() == 0 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 100; ++i) 
        v.append(i);
    REQUIRE( v.getRows() == 100 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 100; ++i)
        REQUIRE( v[i] == i );    
}

TEST_CASE( "Vector resize", "[vector]" ) {
    mp::vector v;
    REQUIRE( v.getRows() == 0 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 100; ++i) 
    {
        v.resize(i);
        REQUIRE( v.getRows() == i );
        REQUIRE( v.getCols() == 1 );
        for (int j = 0; j != i; ++j)
            REQUIRE( v[j] == 0 );
    }
    for (int i = 100; i != -1; --i) 
    {
        v.resize(i);
        REQUIRE( v.getRows() == i );
        REQUIRE( v.getCols() == 1 );
        for (int j = 0; j != i; ++j)
            REQUIRE( v[j] == 0 );
    }  
}

TEST_CASE( "Vector construction - size only", "[vector]" ) {
    mp::vector v(10);
    REQUIRE( v.getRows() == 10 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 10; ++i)
        REQUIRE( v[i] == 0.0 );
}

TEST_CASE( "Vector construction - size and fill value", "[vector]" ) {
    mp::vector v(10, 1.0);
    REQUIRE( v.getRows() == 10 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 10; ++i)
        REQUIRE( v[i] == 1.0 );
}

TEST_CASE( "Vector construction - copy", "[vector]" ) {
    mp::vector w(10, 1.0);
    mp::vector v(w);
    REQUIRE( v.getRows() == 10 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 10; ++i)
        REQUIRE( v[i] == 1.0 );
}

TEST_CASE( "Vector construction - move", "[vector]" ) {
    mp::vector v(mp::vector(10, 1.0));
    REQUIRE( v.getRows() == 10 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 10; ++i)
        REQUIRE( v[i] == 1.0 );    
}

TEST_CASE( "Vector assignment - copy", "[vector]" ) {
    mp::vector w(10, 1.0);
    mp::vector v;
    v = w;
    REQUIRE( v.getRows() == 10 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 10; ++i)
        REQUIRE( v[i] == 1.0 );    
}

TEST_CASE( "Vector assignment - move", "[vector]" ) {
    mp::vector v;
    v = mp::vector(10, 1.0);
    REQUIRE( v.getRows() == 10 );
    REQUIRE( v.getCols() == 1 );
    for (int i = 0; i != 10; ++i)
        REQUIRE( v[i] == 1.0 );    
}

TEST_CASE( "Operations between vectors and scalars/unary operators: +-*/, +-", "[vector]" ) {
    mp::vector v(10);
    mp::vector w = v + 1.0;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 1.0 );    
    w = w - 2.0;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == -1.0 );    
    w = w * (-10.0);
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 10.0 );    
    w = w / 2.0;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 5.0 );
    w = 5.0 + w;   
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 10.0 );    
    w = 5.0 - w;   
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == -5.0 );    
    w = (-1.0) * w;   
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 5.0 );    
    w = 10.0 / w;   
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 2.0 );
    w = -w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == -2.0 );
    w = +w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == -2.0 );
}

TEST_CASE( "Operations between vectors: +-*/", "[vector]" ) {
    mp::vector v(10, 2.0);
    mp::vector w(10, 4.0);
    mp::vector res;

    res = v + w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == 6.0 );
    res = v - w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == -2.0 );
    res = v * w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == 8.0 );
    res = v / w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == .5 );
}

TEST_CASE( "Assignment operators between vectors and scalars: += -= *= /=", "[vector]" ) {
    mp::vector w(10);
    w += 1.0;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 1.0 );    
    w -= 2.0;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == -1.0 );    
    w *= -10.0;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 10.0 );    
    w /= 2.0;
    for (int i = 0; i != 10; ++i)
        REQUIRE( w[i] == 5.0 );    
}

TEST_CASE( "Assignment operators between vectors: += -= *= /=", "[vector]" ) {
    mp::vector v(10, 2.0);
    mp::vector w(10, 4.0);
    mp::vector res;

    res = v;
    res += w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == 6.0 );
    res = v;
    res -= w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == -2.0 );
    res = v;
    res *= w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == 8.0 );
    res = v;
    res /= w;
    for (int i = 0; i != 10; ++i)
        REQUIRE( res[i] == .5 );
}

TEST_CASE( "Vector comparisons: ==, !=" ) {
    mp::vector v(10, 2.0);
    mp::vector w(10, 2.0);
    mp::vector u(10, 1.0);
    mp::vector t(15, 2.0);
    REQUIRE( v == w );
    REQUIRE( v != u );
    REQUIRE( v != t );
}

TEST_CASE( "Vector dot products" ) {
    mp::vector v(10, 1.0);
    mp::vector w(10, 1.0);
    REQUIRE( v.dot(w) == 10.0 );
    REQUIRE( mp::dot(v, w) == 10.0 );
}

TEST_CASE( "Vector math functions" ) {
    mp::vector v(10, -2.0);
    REQUIRE( mp::abs(v) == mp::vector(10, std::abs(-2.0)) );
    REQUIRE( mp::exp(v) == mp::vector(10, std::exp(-2.0)) );
    REQUIRE( mp::exp2(v) == mp::vector(10, std::exp2(-2.0)) );

    v = mp::vector(10, 2.0);
    REQUIRE( mp::log(v) == mp::vector(10, std::log(2.0)) );
    REQUIRE( mp::log10(v) == mp::vector(10, std::log10(2.0)) );
    REQUIRE( mp::log2(v) == mp::vector(10, std::log2(2.0)) );
    REQUIRE( mp::sqrt(v) == mp::vector(10, std::sqrt(2.0)) );
    REQUIRE( mp::sin(v) == mp::vector(10, std::sin(2.0)) );
    REQUIRE( mp::cos(v) == mp::vector(10, std::cos(2.0)) );
    REQUIRE( mp::tan(v) == mp::vector(10, std::tan(2.0)) );
    v = mp::vector(10, .2);
    REQUIRE( mp::asin(v) == mp::vector(10, std::asin(.2)) );
    REQUIRE( mp::acos(v) == mp::vector(10, std::acos(.2)) );
    REQUIRE( mp::atan(v) == mp::vector(10, std::atan(.2)) );

    v = mp::vector(10, 2.0);
    mp::vector w(10, 4.0);
    REQUIRE( mp::pow(v, 2.0) == w );
    REQUIRE( mp::pow(v, v) == w );
    //w = mp::vector(10, exp(-2.0));
/*
    vector abs(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::abs(v[i]); return w; }
    vector exp(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::exp(v[i]); return w; }
    vector exp2(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::exp2(v[i]); return w; }
    vector log(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::log(v[i]); return w; }
    vector log10(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::log10(v[i]); return w; }
    vector log2(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::log2(v[i]); return w; }
    vector pow(const vector& v, const double& d) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::pow(v[i], d); return w; }
    vector pow(const vector& v1, const vector& v2) { check_dimensions(v1, v2); vector w(v1.rows); for (unsigned long i = 0; i != v1.rows; ++i) w[i] = std::pow(v1[i], v2[i]); return w; }
    vector sqrt(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::sqrt(v[i]); return w; }
    vector sin(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::sin(v[i]); return w; }
    vector cos(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::cos(v[i]); return w; }
    vector tan(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::tan(v[i]); return w; }
    vector asin(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::asin(v[i]); return w; }
    vector acos(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::acos(v[i]); return w; }
    vector atan(const vector& v) { vector w(v.rows); for (unsigned long i = 0; i != v.rows; ++i) w[i] = std::atan(v[i]); return w; }
*/
}
