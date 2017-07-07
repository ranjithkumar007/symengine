/**
 *  \file simplify.h
 *  Includes various functions for simplification
 *
 **/

#include <symengine/functions.h>

namespace SymEngine
{
template <typename Func>
RCP<const Basic> bottom_up(const RCP<const Basic> &barg, const Func &f);

RCP<const Basic> trig_simplify_fu(const RCP<const Basic> &arg);
};