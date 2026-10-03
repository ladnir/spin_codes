#pragma once
#include "../packet8_wider24/Packet8Wide24.h"

namespace spin::research::packet8wide24scaling {
using Block=packet8wide24::Block;
using Plan=packet8wide24::Plan;

// Exactly the frozen construction, with independent routing/output store
// policies. Scratch is disjoint and 64-byte aligned, as in the frozen API.
// Streaming calls finish their stores with an sfence before returning.
void reverseRoute(const Block* input,Block* scratch,const Plan&,bool stream);

// Both paths use the exact flat9 multiplier and shared-parity winner.
// A streaming request falls back to cached stores unless output is 64-byte
// aligned. Output may alias the original input, but must not alias scratch.
void outer(const Block* scratch,Block* output,const Plan&,bool stream);
}
