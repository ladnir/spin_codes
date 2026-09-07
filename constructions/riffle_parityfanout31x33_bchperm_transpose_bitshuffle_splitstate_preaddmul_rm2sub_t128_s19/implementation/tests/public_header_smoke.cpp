#include <libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.h>

int main()
{
	osuCrypto::StructuredSpinB256T128S19 spin;
	osuCrypto::StructuredSpinB256T128S19::Workspace workspace;
	return spin.codeBlocks == 2 * spin.messageBlocks &&
		sizeof(workspace) != 0 ? 0 : 1;
}
