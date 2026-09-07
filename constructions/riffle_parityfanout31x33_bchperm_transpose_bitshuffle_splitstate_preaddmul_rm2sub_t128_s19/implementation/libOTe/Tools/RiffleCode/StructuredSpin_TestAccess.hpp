#pragma once

#include "RiffleExactPermFieldCheckpoint.h"
#include "StructuredSpinB256T128S19.h"

namespace osuCrypto
{
	struct StructuredSpinTestAccess
	{
		static const RiffleExactPermFieldCheckpoint& outer(
			const StructuredSpinB256T128S19& spin) noexcept;

		static const RiffleRm2SubS19Transpose& inner(
			const StructuredSpinB256T128S19& spin) noexcept;
	};
}
