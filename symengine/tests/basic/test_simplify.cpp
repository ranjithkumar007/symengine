#include "catch.hpp"
#include <iostream>

#include <symengine/simplify.cpp>
#include <symengine/add.h>

using SymEngine::csc;
using SymEngine::sec;
using SymEngine::sin;
using SymEngine::symbol;
using SymEngine::mul;
using SymEngine::add;
using SymEngine::one;
using SymEngine::integer;
using SymEngine::neg;
using SymEngine::pi;
using SymEngine::zero;
using SymEngine::pow;
using SymEngine::rational;

TEST_CASE("count trigs", "Simplify")
{
    auto x = symbol("x");
    REQUIRE(count_trigs(add(cos(x), sin(x))) == 2);
    REQUIRE(count_trigs(add(sec(x), mul(tan(x), sin(x)))) == 3);
}

TEST_CASE("TRs", "Simplify")
{
    auto x = symbol("x");
    auto i2 = integer(2), i3 = integer(3);

    // TR1
    REQUIRE(eq(*TR1(one), *one));
    REQUIRE(eq(*TR1(add(mul(i2, csc(x)), sec(x))),
               *add(div(i2, sin(x)), div(one, cos(x)))));
    REQUIRE(eq(*TR1(sec(x)), *div(one, cos(x))));

    // TR2
    REQUIRE(eq(*TR2(one), *one));
    REQUIRE(eq(*TR2(tan(x)), *div(sin(x), cos(x))));
    REQUIRE(eq(*TR2(cot(x)), *div(cos(x), sin(x))));

    // TR3
    REQUIRE(eq(*TR3(one), *one));
    REQUIRE(eq(*TR3(cos(add(div(pi, i2), x))), *neg(sin(x))));
    REQUIRE(eq(*TR3(cos(add(mul(pi, integer(15)), x))), *neg(cos(x))));

    // TR4
    REQUIRE(eq(*TR4(one), *one));
    REQUIRE(eq(*TR4(cos(div(pi, i2))), *zero));
    REQUIRE(eq(*TR4(sin(div(pi, i2))), *one));

    // TR5
    REQUIRE(eq(*TR5(one), *one));
    REQUIRE(eq(*TR5(mul(sin(x), sin(x))), *sub(one, mul(cos(x), cos(x)))));
    REQUIRE(eq(*TR5(mul(sin(mul(x, x)), sin(mul(x, x)))),
               *sub(one, mul(cos(mul(x, x)), cos(mul(x, x))))));
    REQUIRE(eq(*TR5(mul(cos(x), cos(x))), *mul(cos(x), cos(x))));
    REQUIRE(eq(*TR6(pow(sin(x), i3)), *pow(sin(x), i3)));

    // TR6
    REQUIRE(eq(*TR6(one), *one));
    REQUIRE(eq(*TR6(mul(sin(x), sin(x))), *mul(sin(x), sin(x))));
    REQUIRE(eq(*TR6(mul(cos(mul(x, x)), cos(mul(x, x)))),
               *sub(one, mul(sin(mul(x, x)), sin(mul(x, x))))));
    REQUIRE(eq(*TR6(mul(cos(x), cos(x))), *sub(one, mul(sin(x), sin(x)))));
    REQUIRE(eq(*TR6(pow(cos(x), integer(4))),
               *pow(sub(one, mul(sin(x), sin(x))), i2)));
    REQUIRE(eq(*TR6(pow(cos(x), i3)), *pow(cos(x), i3)));

    // TR7
    REQUIRE(eq(*TR7(one), *one));
    REQUIRE(eq(*TR7(mul(cos(x), cos(x))), *div(add(one, cos(mul(i2, x))), i2)));
    REQUIRE(eq(*TR7(add(cos(x), mul(cos(x), cos(x)))),
               *add(cos(x), div(add(one, cos(mul(i2, x))), i2))));

    // TR8
    REQUIRE(eq(*TR8(one), *one));
    REQUIRE(eq(*TR8(mul(cos(i2), cos(i3))),
               *add(div(cos(one), i2), div(cos(integer(5)), i2))));
    REQUIRE(eq(*TR8(mul(cos(i2), sin(i3))),
               *add(div(sin(one), i2), div(sin(integer(5)), i2))));
    REQUIRE(eq(*TR8(mul(sin(i2), sin(i3))),
               *sub(div(cos(one), i2), div(cos(integer(5)), i2))));
    REQUIRE(eq(*TR8(mul({sin(one), sin(i2), sin(i3)})),
               *add({div(sin(i2), integer(4)), div(sin(integer(4)), integer(4)),
                     div(sin(integer(6)), integer(-4))})));
    REQUIRE(eq(*TR8(mul({cos(i2), cos(i3), cos(integer(4)), cos(integer(5))})),
               *add({div(one, integer(8)), div(cos(i2), integer(8)),
                     div(cos(integer(4)), integer(4)),
                     div(cos(integer(6)), integer(8)),
                     div(cos(integer(8)), integer(8)),
                     div(cos(integer(10)), integer(8)),
                     div(cos(integer(14)), integer(8))})));
    REQUIRE(eq(*TR8(mul({cos(i2), cos(i3), cos(integer(4)), cos(integer(5)),
                         cos(integer(6))})),
               *add({div(one, integer(16)), mul(cos(i2), rational(3, 16)),
                     div(cos(integer(4)), integer(8)),
                     div(cos(integer(6)), integer(8)),
                     div(cos(integer(8)), integer(8)),
                     div(cos(integer(10)), integer(8)),
                     div(cos(integer(12)), integer(16)),
                     div(cos(integer(14)), integer(16)),
                     div(cos(integer(16)), integer(16)),
                     div(cos(integer(20)), integer(16))})));
    REQUIRE(eq(*TR8(mul({pow(sin(div(mul(i3, pi), integer(7))), i2),
                         pow(cos(div(mul(i3, pi), integer(7))), i2)})),
               *sub(div(one, integer(8)),
                    div(sin(mul(rational(3, 14), pi)), integer(8)))));


    // TR9
    auto half = div(one, i2);
    REQUIRE(eq(*TR9(one), *one));
    REQUIRE(eq(*TR9(add(sin(x), cos(x))), *add(sin(x), cos(x))));
    REQUIRE(eq(*TR9(add(cos(one), cos(i2))),
               *mul({i2, cos(half), cos(mul(i3, half))})));
    REQUIRE(eq(*TR9(sub(cos(one), cos(i2))),
               *mul({i2, sin(half), sin(mul(i3, half))})));
    REQUIRE(eq(*TR9(sub(sin(one), sin(i2))),
               *neg(mul({i2, sin(half), cos(mul(i3, half))}))));
    REQUIRE(eq(*TR9(add(sin(one), sin(i2))),
               *mul({i2, cos(half), sin(mul(i3, half))})));
    REQUIRE(
        eq(*TR9(add({cos(i2), cos(integer(4)), mul({i2, cos(one), cos(i3)})})),
           *mul({integer(4), cos(one), cos(i3)})));
}

TEST_CASE("fu", "Simplify")
{
    auto x = symbol("x");
    auto i2 = integer(2);

    std::cout << *trig_simplify_fu(sec(x)) << "\n";
    std::cout << *trig_simplify_fu(add(sec(x), csc(x))) << "\n";
    std::cout << *trig_simplify_fu(add(sec(mul(x, i2)), mul(sec(x), sin(x))))
              << "\n";
    std::cout << *trig_simplify_fu(
                     add(mul(cos(x), cos(x)), mul(sin(x), sin(x))))
              << "\n";
}
