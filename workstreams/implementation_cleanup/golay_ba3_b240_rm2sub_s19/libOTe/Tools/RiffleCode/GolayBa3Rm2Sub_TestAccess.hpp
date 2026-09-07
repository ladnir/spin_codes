#pragma once

#include <libOTe/Tools/RiffleCode/GolayBa3Rm2SubB240K20.h>

namespace osuCrypto
{
	struct GolayBa3Rm2SubTestAccess
	{
		static void validateSetup(const GolayBa3Rm2SubB240K20& code);

		static void dualEncodeReference(
			const GolayBa3Rm2SubB240K20& code,
			const block* input,
			block* output,
			GolayBa3Rm2SubB240K20::Workspace& workspace);
	};
}
