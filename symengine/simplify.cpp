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
    if (is_a<Add>(*barg) or is_a<Mul>(*barg)) {
        vec_basic newargs;
        for (const auto &a : barg->get_args()) {
            newargs.push_back(bottom_up(a, f));
        }
        if (is_a<Add>(*barg))
            return f(add(newargs));
        return f(mul(newargs));
    }
    if (is_a<Pow>(*barg)) {
        auto &dfunc = down_cast<const Pow &>(*barg);
        auto farg1 = dfunc.get_base(), farg2 = dfunc.get_exp();
        auto newarg1 = bottom_up(farg1, f), newarg2 = bottom_up(farg2, f);
        auto nbarg = barg;
        if (farg1 != newarg1 or farg2 != newarg2)
            nbarg = pow(newarg1, newarg2);
        return f(nbarg);
    }
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
    int num_ops = 0;
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

    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic> {
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

    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic> {
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

    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic> {
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
                return pow(sub(one, pow(cos(b->get_args()[0]), i2)),
                           div(e, i2));
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

    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic> {
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
                return pow(sub(one, pow(sin(b->get_args()[0]), i2)),
                           div(e, i2));
            }
        }
        return a;
    };

    return bottom_up(arg, f);
}

RCP<const Basic> TR7(const RCP<const Basic> &arg)
{
    // Lowering the degree of cos’ square.

    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic> {
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

    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic> {
        if (not(is_a<Pow>(*a)
                and is_a<Integer>(*down_cast<const Pow &>(*a).get_exp())
                and (is_a<Sin>(*down_cast<const Pow &>(*a).get_base())
                     or is_a<Cos>(*down_cast<const Pow &>(*a).get_base())))
            and not is_a<Mul>(*a))
            return a;
        vec_basic cos_args, sin_args, others;
        for (const auto &e : a->get_args()) {
            if (is_a<Cos>(*e)) {
                cos_args.push_back(e->get_args()[0]);
            } else if (is_a<Sin>(*e)) {
                sin_args.push_back(e->get_args()[0]);
            } else if (is_a<Pow>(*e)) {
                if (is_a<Sin>(*down_cast<const Pow &>(*e).get_base())
                    and is_a<Integer>(*down_cast<const Pow &>(*e).get_exp())
                    and eq(*Gt(down_cast<const Pow &>(*e).get_exp(), zero),
                           *boolTrue)) {
                    auto n = down_cast<const Pow &>(*e).get_exp();
                    while (eq(*Gt(n, zero), *boolTrue)) {
                        sin_args.push_back(down_cast<const Pow &>(*e)
                                               .get_base()
                                               ->get_args()[0]);
                        n = sub(n, one);
                    }
                } else if (is_a<Cos>(*down_cast<const Pow &>(*e).get_base())
                           and is_a<Integer>(
                                   *down_cast<const Pow &>(*e).get_exp())
                           and eq(*Gt(down_cast<const Pow &>(*e).get_exp(),
                                      zero),
                                  *boolTrue)) {
                    auto n = down_cast<const Pow &>(*e).get_exp();
                    while (eq(*Gt(n, zero), *boolTrue)) {
                        cos_args.push_back(down_cast<const Pow &>(*e)
                                               .get_base()
                                               ->get_args()[0]);
                        n = sub(n, one);
                    }
                } else {
                    others.push_back(e);
                }
            } else {
                others.push_back(e);
            }
        }
        if ((sin_args.empty() and cos_args.size() <= 1)
            or (cos_args.empty() and sin_args.size() <= 1))
            return a;
        auto n = std::min(sin_args.size(), cos_args.size());
        RCP<const Basic> alpha, beta;
        auto i2 = integer(2);

        while (n--) {
            // sin α * cos β =  1/2 [sin(α + β) + sin(α − β)]
            alpha = sin_args.back(), beta = cos_args.back();
            sin_args.pop_back();
            cos_args.pop_back();
            others.push_back(
                div(add(sin(add(alpha, beta)), sin(sub(alpha, beta))), i2));
        }
        while (cos_args.size() > 1) {
            // cos α · cos β = 1/2 [cos(α + β) + cos(α − β)]
            alpha = cos_args.back();
            cos_args.pop_back();
            beta = cos_args.back();
            cos_args.pop_back();
            others.push_back(
                div(add(cos(add(alpha, beta)), cos(sub(alpha, beta))), i2));
        }
        if (not cos_args.empty())
            others.push_back(cos(cos_args.back()));
        while (sin_args.size() > 1) {
            // sin α · sin β = − 1/2 [cos(α + β) − cos(α − β)]
            alpha = sin_args.back();
            sin_args.pop_back();
            beta = sin_args.back();
            sin_args.pop_back();
            others.push_back(div(
                sub(cos(add(alpha, beta)), cos(sub(alpha, beta))), neg(i2)));
        }
        if (not sin_args.empty())
            others.push_back(sin(sin_args.back()));
        return TR8(expand(mul(others)));
    };

    return bottom_up(arg, f);
}

RCP<const Basic> TR9(const RCP<const Basic> &arg)
{
    // Converting sum or difference to product:
    // eg : sin α + sin β = 2 * sin (α+β)/2 * cos (α-β)/2
    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic> {
        if (not(is_a<Add>(*a) and a->get_args().size() > 1))
            return a;
        vec_basic sin_args, cos_args, others;
        // TODO : presently, this can't handle inputs of form a*sin(x)(or
        // cos(x)) + b*sin(y)(or cos(y)) where a != +/- 1,b != +/-1.
        for (const auto &e : a->get_args()) {
            if (is_a<Sin>(*e)) {
                sin_args.push_back(e->get_args()[0]);
            } else if (is_a<Cos>(*e)) {
                cos_args.push_back(e->get_args()[0]);
            } else if (is_a<Mul>(*e)
                       and eq(*down_cast<const Mul &>(*e).get_coef(),
                              *minus_one)) {
                auto negE = neg(e);
                if (is_a<Sin>(*negE)) {
                    sin_args.push_back(neg(negE->get_args()[0]));
                } else if (is_a<Cos>(*negE)) {
                    cos_args.push_back(add(pi, negE->get_args()[0]));
                } else {
                    others.push_back(e);
                }
            } else {
                others.push_back(e);
            }
        }
        if (sin_args.size() <= 1 and cos_args.size() <= 1)
            return a;
        RCP<const Basic> alpha, beta;
        auto i2 = integer(2);
        if (sin_args.size() > 1) {
            for (size_t i = 0; i < sin_args.size() - 1; i += 2) {
                // sin α + sin β = 2 * sin (α+β)/2 * cos (α-β)/2
                others.push_back(mul(
                    {i2,
                     sin(expand(div(add(sin_args[i], sin_args[i + 1]), i2))),
                     cos(expand(div(sub(sin_args[i], sin_args[i + 1]), i2)))}));
            }
        }
        if (cos_args.size() > 1) {
            for (size_t i = 0; i < cos_args.size() - 1; i += 2) {
                // cos α + cos β = 2 * cos (α+β)/2 · cos (α−β)/2
                others.push_back(expand(mul(
                    {i2,
                     cos(expand(div(add(cos_args[i], cos_args[i + 1]), i2))),
                     cos(expand(
                         div(sub(cos_args[i], cos_args[i + 1]), i2)))})));
            }
        }
        return TR9(expand(add(others)));
    };

    return bottom_up(arg, f);
}

RCP<const Basic> TR10(const RCP<const Basic> &arg)
{
    // Sum or difference of angles
    // eg : sin(α + β) = sin α cos β + cos α sin β
    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic>
    {
        if(not (is_a<Sin>(*a) or is_a<Cos>(*a)))
            return a;
        auto arg = a->get_args()[0]
        if (not is_a<Add>(*arg))
            return a;
        auto addArgs = arg->get_args();
        auto a = addArgs.back();
        addArgs.pop_back();
        if (addArgs.size() == 1) {
            auto b = addArgs[0];
            if (is_a<Sin>(*a)) {
                return add(mul(sin(a), cos(b)), mul(cos(a), sin(b)));
            } else {
                return sub(mul(cos(a), cos(b)), mul(sin(a), sin(b)));
            }
        } else {
            auto b = add(addArgs);
            if (is_a<Sin>(*a)) {
                return add(mul(sin(a), TR10(cos(b))), mul(cos(a), TR10(sin(b))));
            } else {
                return sub(mul(cos(a), TR10(cos(b))), mul(sin(a), TR10(sin(b))));
            }
        }
    }

    return bottom_up(arg, f);
}

RCP<const Basic> TR11(const RCP<const Basic> &arg)
{
    // Double angle formulas:
    // eg : sin 2α = 2 sin α cos α
    static auto f = [](const RCP<const Basic> &a) -> RCP<const Basic>
    {
        if(not (is_a<Sin>(*a) or is_a<Cos>(*a)))
            return a;
        if (not is_a_Number(*a)) {
            auto arg = div(a.get_args()[0], integer(2));
            auto cval = cos(arg), sval = sin(arg);
            if (is_a<Sin>(*a)) {
                // sin 2α = 2 sin α cos α
                return mul({i2, cval, sval});
            } else {
                // cos 2α = (cos(α))**2 − (sin(α))**2
                return sub(pow(cval, i2), pow(sval, i2));
            }
        }
        return a;
    }

    return bottom_up(arg, f);
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
        = {std::numeric_limits<int>::max(), zero}; // init to some dummy node.
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
    auto narg = TR1(arg);
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
