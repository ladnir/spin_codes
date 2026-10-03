#pragma once
#include "../packet8_wider24/Packet8Wide24.h"
namespace spin::research::packet8wide24opt {
namespace base=spin::research::packet8wide24;
using base::Block;
using base::Plan;
void reverseRouteTraffic(const Block*,Block*,const Plan&,unsigned variant);
}
