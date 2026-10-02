#pragma once
#include "RsWide.h"
namespace spin::research::rswide {
void outerFastPair(const Block*,Block*,const Plan&);
void outerFastNoInline(const Block*,Block*,const Plan&);
void outerFastPairNoInline(const Block*,Block*,const Plan&);
void transposeFastPair(const Block*,Block*,Block*,const Plan&);
void transposeFastNoInline(const Block*,Block*,Block*,const Plan&);
void transposeFastPairNoInline(const Block*,Block*,Block*,const Plan&);
void outerFastPlaneNoInline(const Block*,Block*,const Plan&);
void outerFastFusedPair(const Block*,Block*,const Plan&);
void transposeFastPlaneNoInline(const Block*,Block*,Block*,const Plan&);
void transposeFastFusedPair(const Block*,Block*,Block*,const Plan&);
void outerFastFusedPairNt(const Block*,Block*,const Plan&);
void transposeFastFusedPairNt(const Block*,Block*,Block*,const Plan&);
}
