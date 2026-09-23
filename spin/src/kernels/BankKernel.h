#pragma once
#include <spin/detail/BankState.h>
#include "Block.h"
namespace spin::detail::kernel {
void bankTranspose(const BankState&,bool four,const void* in,void* out,storage::block* scratch,std::uint32_t* addresses);
// False means a rejected random word: regenerate the complete scalar stream.
bool bankMasks512(std::uint32_t* masks,std::size_t count,unsigned bits,std::uint64_t seed) noexcept;
}
