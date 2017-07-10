#include <symengine/simplify.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/logic.h>

namespace SymEngine
{
template <typename Func>
RCP<const Basic> bottom_up(const RCP<const Basic> &barg, const Func &f)
{
    if (is_a_sub<OneArgFunction>(*barg)) {
        const OneArgFunction &dfunc = down_cast<const OneArgFunction &>(*barg);
        auto farg = dfunc.get_arg();
        auto newarg = bottom_up(farg, f);
        auto nbarg = barg;
        if (not eq(*newarg, *farg))
            nbarg = dfunc.create(newarg);
        return f(nbarg);
    }
    if (is_a_sub<TwoArgFunction>(*barg)) {
        const TwoArgFunction &dfunc = down_cast<const TwoArgFunction &>(*barg);
        auto farg1 = dfunc.get_arg1(), farg2 = dfunc.get_arg2();
        auto newarg1 = bottom_up(farg1, f), newarg2 = bottom_up(farg2, f);
        auto nbarg = barg;
        if (farg1 != newarg1 or farg2 != newarg2)
            nbarg = dfunc.create(newarg1, newarg2);
        return f(nbarg);
    }
    if (is_a_sub<MultiArgFunction>(*barg)) {
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
    if (is_a_sub<TrigFunction>(*arg))
        return (1 + count_trigs(arg->get_args()[0]));

    for (const auto &a : arg->get_args()) {
        num_ops += count_trigs(a);
    }
    return num_ops;
}

inline RCP<const Basic> identity(const RCP<const Basic> &arg)
{
    return arg;
}

inline RCP<const Basic> TR0(const RCP<const Basic> &arg)
{
    // Simplify the expression
    return expand(arg);
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

inline RCP<const Basic> TR3(const RCP<const Basic> &arg)
{
    // Induced formulas like sin(2*k*pi + x) = sin(x)
    return arg; // automatically handled. for all trig functions.
}

inline RCP<const Basic> TR4(const RCP<const Basic> &arg)
{
    // Values of special angles (0, pi/6, pi/4, pi/3, pi/2)
    return arg; // automatically handled.
}

RCP<const Basic> TR5(const RCP<const Basic> &arg)
{
    // Substitution of sin’s square:
    // (sin(x))**2 = 1 − (cos(x))**2.

    auto f = [](const RCP<const Basic> &a) {
        if (is_a<Pow>(*a)
            and is_a<Sin>(*down_cast<const Pow &>(*a).get_base())) {
            auto b = down_cast<const Pow &>(*a).get_base();
            auto e = down_cast<const Pow &>(*a).get_exp();
            if (eq(*Ge(e, zero), *boolFalse))
                return a;
            if (is_a<Integer>(*e)) {
                // currently mod() supports only integers.
                auto i2 = integer(2);
                if (not eq(*mod(down_cast<const Integer &>(*e), *i2), *zero))
                    return a;
                return pow(sub(one, pow(cos(b), i2)), div(e, i2));
            }
        }
        return a;
    };

    return bottom_up(arg, f);
}

RCP<const Basic> TR6(const RCP<const Basic> &arg)
{
    // Substitution of cos's square:
    // (cos(x))**2 = 1 − (sin(x))**2.

    auto f = [](const RCP<const Basic> &a) {
        if (is_a<Pow>(*a)
            and is_a<Cos>(*down_cast<const Pow &>(*a).get_base())) {
            auto b = down_cast<const Pow &>(*a).get_base();
            auto e = down_cast<const Pow &>(*a).get_exp();
            if (eq(*Ge(e, zero), *boolFalse))
                return a;
            if (is_a<Integer>(*e)) {
                // currently mod() supports only integers.
                auto i2 = integer(2);
                if (not eq(*mod(down_cast<const Integer &>(*e), *i2), *zero))
                    return a;
                return pow(sub(one, pow(sin(b), i2)), div(e, i2));
            }
        }
        return a;
    };

    return bottom_up(arg, f);
}

RCP<const Basic> TR7(const RCP<const Basic> &arg)
{
    // Lowering the degree of cos’ square.
    auto f = [](const RCP<const Basic> &a) {
        auto i2 = integer(2);
        if (is_a<Pow>(*a) and is_a<Cos>(*down_cast<const Pow &>(*a).get_base())
            and eq(*down_cast<const Pow &>(*a).get_exp(), *i2)) {
            return div(
                add(one, cos(mul(i2, down_cast<const Cos &>(
                                         *down_cast<const Pow &>(*a).get_base())
                                         .get_arg()))),
                i2);
        }
        return a;
    };
    return bottom_up(arg, f);
}

RCP<const Basic> TR8(const RCP<const Basic> &arg)
{
    // Converting product to sum or difference
    // eg : sin α * cos β =  1/2 [sin(α + β) + sin(α − β)]
    return arg;
}

RCP<const Basic> TR9(const RCP<const Basic> &arg)
{
    // Converting sum or difference to product:
    // eg : sin α + sin β = 2 sin (α+β)/2 * cos (α-β)/2
    return arg;
}

RCP<const Basic> TR10(const RCP<const Basic> &arg)
{
    // Sum or difference of angles
    // eg : sin(α + β) = sin α cos β + cos α sin β
    return arg;
}

RCP<const Basic> TR11(const RCP<const Basic> &arg)
{
    // Double angle formulas:
    // eg : sin 2α = 2 sin α cos α
    return arg;
}

RCP<const Basic> TR12(const RCP<const Basic> &arg)
{
    // Sum or difference of tan:
    // eg : tan(α + β) = (tan α+tan β) / (1− (tanα * tan β))
    return arg;
}

RCP<const Basic> TR13(const RCP<const Basic> &arg)
{
    // Product of tan or cot:
    // eg : tan α * tan β = 1 − (tan α + tan β) * cot(α + β)
    return arg;
}

inline std::pair<int, RCP<const Basic>>
compare_pair(const std::pair<int, RCP<const Basic>> &a,
             const std::pair<int, RCP<const Basic>> &b)
{
    if (a.first <= b.first) {
        return a;
    }
    return b;
}

typedef std::vector<std::vector<std::function<RCP<const Basic>(
    const RCP<const Basic> &)>>> vec2D_func_basic_basic;

RCP<const Basic> find_best_choice(const RCP<const Basic> &arg,
                                  const vec2D_func_basic_basic &a)
{
    RCP<const Basic> narg;
    std::pair<int, RCP<const Basic>> bestCh
        = {INT_MAX, zero}; // init to some dummy node.
    for (const auto &choice : a) {
        narg = arg;
        for (const auto &chain : choice) {
            narg = chain(narg);
        }
        bestCh = compare_pair(bestCh, {count_trigs(narg), narg});
    }
    return bestCh.second;
}

RCP<const Basic> CTR1(const RCP<const Basic> &arg)
{
    // CTR1: Combination rule 1
    // Intuitive rules = {TR5, TR6, TR0}
    vec2D_func_basic_basic v = {{identity}, {TR5, TR0}, {TR6, TR0}};
    return find_best_choice(arg, v);
}

RCP<const Basic> CTR2(const RCP<const Basic> &arg)
{
    // CTR2: Combination rule 2.
    // Intuitive rules = {TR5, TR6, TR11}
    vec2D_func_basic_basic v = {{TR0}, {TR5, TR0}, {TR6, TR0}};
    auto ch = TR11(arg);
    return find_best_choice(ch, v);
}

RCP<const Basic> CTR3(const RCP<const Basic> &arg)
{
    // CTR3: Combination rule 3.
    // Intuitive rules = {TR8, TR10inverse}
    // TODO : {TR8, TR10inv, TR0} once TR10inv gets implemented.

    vec2D_func_basic_basic v = {{TR8, TR0}, {identity}};
    return find_best_choice(arg, v);
}

RCP<const Basic> CTR4(const RCP<const Basic> &arg)
{
    // CTR4: Combination rule 4.
    // Intuitive rules = {TR4, TR10inverse}
    // TODO : {TR4, TR10inv} once TR10inv gets implemented.
    vec2D_func_basic_basic v = {{identity}};
    return find_best_choice(arg, v);
}

RCP<const Basic> apply_RL1(const RCP<const Basic> &arg)
{
    // RL1: Rule List 1
    vec2D_func_basic_basic v
        = {{identity, TR4, TR3, TR4, TR12, TR4, TR13, TR4, TR0}};
    return find_best_choice(arg, v);
}

RCP<const Basic> apply_RL2(const RCP<const Basic> &arg)
{
    // RL2: Rule List 2
    vec2D_func_basic_basic v
        = {{TR4, TR3, TR10, TR4, TR3, TR11},
           {TR5, TR7, TR11, TR4},
           {CTR3, CTR1, TR9, CTR2, TR4, TR9, TR0, TR9, CTR4},
           {identity}};
    return find_best_choice(arg, v);
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

    auto narg = std::move(TR1(arg));
    if (has<Tan>(narg) or has<Cot>(narg)) {
        auto targ = apply_RL1(narg);
        if (count_trigs(targ) < count_trigs(narg))
            narg = targ;
    }

    if (has<Tan>(narg) or has<Cot>(narg)) {
        narg = TR2(narg);
    }

    narg = TR0(narg);

    if (has<Sin>(narg) or has<Cos>(narg)) {
        auto targ = apply_RL2(narg);
        if (count_trigs(targ) < count_trigs(narg))
            narg = targ;
    }
    return narg;
}
}
