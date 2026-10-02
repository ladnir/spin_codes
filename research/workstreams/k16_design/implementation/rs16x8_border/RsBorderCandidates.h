#pragma once
#include "RsBorder.h"
namespace spin::research::rsborder {
void reverseRouteFused(const Block*,Block*,const Plan&);
void transposeFastFused(const Block*,Block*,Block*,const Plan&);
void reverseRouteTernary(const Block*,Block*,const Plan&);
void transposeFastTernary(const Block*,Block*,Block*,const Plan&);
void reverseRouteEarly(const Block*,Block*,const Plan&);
void transposeFastEarly(const Block*,Block*,Block*,const Plan&);
}
