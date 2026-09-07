#include <libOTe/Tools/RiffleCode/GolayBa3Rm2SubB240K20.h>

int main()
{
	osuCrypto::GolayBa3Rm2SubB240K20 code;
	osuCrypto::GolayBa3Rm2SubB240K20::Workspace workspace;
	return code.codeBlocks == 2'119'680 &&
		code.messageBlocks == (osuCrypto::u64{ 1 } << 20) &&
		code.shortenedMessageBlocks == 11'264 &&
		code.innerEpochs == 16'560 ? 0 : 1;
}
