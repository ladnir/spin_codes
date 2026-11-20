#pragma once
#include <spin/Code.h>
#include "rs8/Rs8.h"

namespace spin::detail {
// Resolves a requested backend for this CPU and build; throws if unavailable.
Backend selectTurboBackend(Backend requested);
// Setup path that produces the selected backend's tables (all bit-identical).
certified22::SetupPath turboSetupPath(Backend) noexcept;
// Transposed encoding 2K -> K with the selected kernel. Scratch holds
// plan.scratchBlocks() blocks; exact input/output aliasing is allowed.
void turboTranspose(Backend, const storage::block* input, storage::block* output,
                    storage::block* scratch, const rs8::Plan& plan);
}
