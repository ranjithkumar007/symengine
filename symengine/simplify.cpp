#include <symengine/simplify.h>
#include <symengine/add.h>
#include <symengine/mul.h>

namespace SymEngine
{
template <typename Func>
RCP<const Basic> bottom_up(const RCP<const Basic> &barg, const Func &f)
{
    if (is_a_OneArgFunction(*barg)) {
        const OneArgFunction &dfunc = down_cast<const OneArgFunction &>(*barg);
        auto farg = dfunc.get_arg();
        auto newarg = bottom_up(farg, f);
        auto nbarg = barg;
        if (not eq(*newarg, *farg))
            nbarg = dfunc.create(newarg);
        return f(nbarg);
    }
    if (is_a<ATan2>(*barg) or is_a<LowerGamma>(*barg) or is_a<UpperGamma>(*barg)
        or is_a<Beta>(*barg)) {
        const TwoArgFunction &dfunc = down_cast<const TwoArgFunction &>(*barg);
        auto farg1 = dfunc.get_arg1(), farg2 = dfunc.get_arg2();
        auto newarg1 = bottom_up(farg1, f), newarg2 = bottom_up(farg2, f);
        auto nbarg = barg;
        if (farg1 != newarg1 or farg2 != newarg2)
            nbarg = dfunc.create(newarg1, newarg2);
        return f(nbarg);
    }
    if (is_a<LeviCivita>(*barg) or is_a<Max>(*barg) or is_a<Min>(*barg)) {
        const MultiArgFunction &dfunc
            = down_cast<const MultiArgFunction &>(*barg);
        auto fargs = dfunc.get_args();
        vec_basic newargs;
        for (const auto &a : fargs) {
            newargs.push_back(bottom_up(a, f));
        }
        auto nbarg = dfunc.create(newargs);
        return f(nbarg);
    }
    return barg;
}

int count_trigs(const RCP<const Basic> &arg)
{
    // acts as measure for choosing a `simplified` form.
    int num_ops = 1;
    if (is_a_TrigFunction(*arg))
        return (1 + count_trigs(arg->get_args()[0]));

    for (const auto &a : arg->get_args()) {
        num_ops += count_trigs(a);
    }
    return num_ops;
}

RCP<const Basic> TR1(const RCP<const Basic> &arg)
{
    // Replaces sec with 1/cos and csc with 1/sin

    auto f = [](const RCP<const Basic> &a) {
        if (is_a<Sec>(*a)) {
            return div(one, cos(a->get_args()[0]));
        }
        if (is_a<Csc>(*a)) {
            return div(one, sin(a->get_args()[0]));
        }
        return a;
    };
    return bottom_up(arg, f);
}

RCP<const Basic> TR2(const RCP<const Basic> &arg)
{
    // Replaces tan with sin/cos and cot with cos/sin

    auto f = [](const RCP<const Basic> &a) {
        if (is_a<Tan>(*a)) {
            return div(sin(a->get_args()[0]), cos(a->get_args()[0]));
        }
        if (is_a<Cot>(*a)) {
            return div(cos(a->get_args()[0]), sin(a->get_args()[0]));
        }
        return a;
    };
    return bottom_up(arg, f);
}

RCP<const Basic> apply_RL1(const RCP<const Basic> &arg)
{
    // TR4 -> TR3 -> TR4 -> TR12 -> TR4 -> TR13 -> TR4 -> TR0
    return arg;
}

RCP<const Basic> apply_RL2(const RCP<const Basic> &arg)
{
    // TR4 -> TR3 -> TR10 -> TR4 -> TR3 -> TR11 -> TR5 -> TR7 -> TR11 -> TR4 ->
    // CTR3 -> CTR1 -> TR9 -> CTR2 -> TR4 -> TR9 -> TR9 -> CTR4
    return arg;
}

template <typename T>
bool has(const RCP<const Basic> &b)
{
    for (const auto &a : b->get_args()) {
        if (is_a<T>(*a))
            return true;
        if (has<T>(a))
            return true;
    }
    return false;
}

// ref for fu :
// https://ai2-s2-pdfs.s3.amazonaws.com/718d/67a8d1ce0a23808c1cc265612f81977357e4.pdf
// See `Figure 4 in the above paper` for Flowchart of the simplification
// procedure.
RCP<const Basic> trig_simplify_fu(const RCP<const Basic> &arg)
{
    if (is_a<Add>(*arg) or is_a<Mul>(*arg)) {
        vec_basic newargs;
        for (const auto &a : arg->get_args())
            newargs.push_back(trig_simplify_fu(a));
        if (is_a<Add>(*arg))
            return add(newargs);
        return mul(newargs);
    }

    // TODO : make a HasVisitor for any general function.
    auto narg = std::move(TR1(arg));
    if (has<Tan>(narg) or has<Cot>(narg)) {
        auto targ = apply_RL1(narg);
        if (count_trigs(targ) < count_trigs(narg))
            narg = targ;
    }

    if (has<Tan>(narg) or has<Cot>(narg)) {
        narg = TR2(narg);
    }

    // narg = TR0(narg);

    if (has<Sin>(narg) or has<Cos>(narg)) {
        auto targ = apply_RL2(narg);
        if (count_trigs(targ) < count_trigs(narg))
            narg = targ;
    }
    return narg;
}
}
