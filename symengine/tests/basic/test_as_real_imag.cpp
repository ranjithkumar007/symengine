#include "catch.hpp"

#include <symengine/real_imag.cpp>

using SymEngine::Basic;
using SymEngine::Add;
using SymEngine::Mul;
using SymEngine::Complex;
using SymEngine::Symbol;
using SymEngine::I;
using SymEngine::sqrt;
using SymEngine::RCP;
using SymEngine::zero;
using SymEngine::neg;
using SymEngine::integer;

TEST_CASE("RealImag: Pow", "[as_real_imag]")
{
    auto sq = sqrt(neg(I));
    RCP<const Basic> re, im;
    auto i2 = integer(2);

    as_real_imag(sqrt(neg(I)), outArg(re), outArg(im));
    REQUIRE(eq(*re, *div(sqrt(i2), i2)));
    REQUIRE(eq(*im, *neg(div(sqrt(i2), i2))));

    as_real_imag(neg(sqrt(neg(I))), outArg(re), outArg(im));
    REQUIRE(eq(*re, *neg(div(sqrt(i2), i2))));
    REQUIRE(eq(*im, *div(sqrt(i2), i2)));
}