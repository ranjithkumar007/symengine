/**
 *  \file simplify.h
 *  Includes various functions for simplification
 *
 **/

#include <symengine/functions.h>
#include <functional>

namespace SymEngine
{
RCP<const Basic> bottom_up(const RCP<const Basic> &rv,
                           const std::function<RCP<const Basic>> &f)
{
    if (is_a_OneArgFunction(*rv)) {
        const OneArgFunction &dfunc = down_cast<const OneArgFunction &>(*rv);
        auto farg = dfunc.get_arg();
        auto newarg = bottom_up(farg, f);
        auto nrv = rv;
        if (not eq(*newarg, *farg))
            nrv = dfunc.create(newarg);
        return f(nrv);
    }
    if (is_a<ATan2>(*rv) or is_a<LowerGamma>(*rv) or is_a<UpperGamma>(*rv)
        or is_a<Beta>(*rv)) {
        const TwoArgFunction &dfunc = down_cast<const TwoArgFunction &>(*rv);
        auto farg1 = dfunc.get_arg1(), farg2 = dfunc.get_arg2();
        auto newarg1 = bottom_up(farg1, f), newarg2 = bottom_up(farg2, f);
        auto nrv = rv;
        if (farg1 != newarg1 or farg2 != newarg2)
            nrv = dfunc.create(newarg1, newarg2);
        return f(nrv);
    }
    if (is_a<LeviCivita>(*rv) or is_a<Max>(*rv) or is_a<Min>(*rv)) {
        const MultiArgFunction &dfunc
            = down_cast<const MultiArgFunction &>(*rv);
        auto fargs = dfunc.get_args();
        vec_basic newargs;
        for (const auto &a : fargs) {
            newargs.push_back(bottom_up(a, f));
        }
        auto nrv = rv;
        if (not eq(*fargs, *newargs)) {
            nrv = dfunc.create(nrv);
        }
        return f(nrv);
    }
    return rv;
}
};