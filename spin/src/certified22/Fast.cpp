#include "Plan.h"

namespace spin::detail::certified22 {
void transposeFast(const Block* input,Block* output,Block* scratch,const Plan& plan) {
    reverseRouteFast(input,scratch,plan);
    if(plan.k<=262144) outerFastCached(scratch,output,plan);
    else outerFastStreaming(scratch,output,plan);
}
}
