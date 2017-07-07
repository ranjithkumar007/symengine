#include "catch.hpp"
#include <iostream>

#include <symengine/simplify.h>
#include <symengine/add.h>

using SymEngine::csc;
using SymEngine::sec;
using SymEngine::sin;
using SymEngine::symbol;
using SymEngine::mul;
using SymEngine::add;
using SymEngine::integer;

TEST_CASE("fu", "Simplify")
{
    auto x = symbol("x");
    auto i2 = integer(2);

    std::cout << *trig_simplify_fu(sec(x)) << "\n";
    std::cout << *trig_simplify_fu(add(sec(x), csc(x))) << "\n";
    std::cout << *trig_simplify_fu(add(sec(mul(x, i2)), mul(sec(x), sin(x))))
              << "\n";
}
