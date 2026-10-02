"""Owned-buffer geometry screen; no encoder or setup distribution changes."""
import argparse
from pathlib import Path
from further_control_codegen import aligned_driver,replace_once


def generate(path):
    source=aligned_driver(path)
    source=replace_once(source,
        'std::vector<k::block> ownedInput(n+12),expected(n),initial(n),guarded((n/1024)*tileStride+12);',
        '''static_assert(SPIN_MEMORY_ALIGNMENT==64 || SPIN_MEMORY_ALIGNMENT==2097152);
    constexpr std::size_t memoryAlignment=SPIN_MEMORY_ALIGNMENT;
    // Both modes reserve equal owner sizes. Only their interior alignment changes.
    constexpr std::size_t slack=2097152/sizeof(k::block)+16;
    std::vector<k::block> ownedInput(n+slack),expected(n),initial(n),guarded((n/1024)*tileStride+slack);''')
    source=replace_once(source,
        '((reinterpret_cast<std::uintptr_t>(ownedInput.data()+4)+63)&~std::uintptr_t(63))+SPIN_INPUT_OFFSET',
        '((reinterpret_cast<std::uintptr_t>(ownedInput.data()+4)+memoryAlignment-1)&~std::uintptr_t(memoryAlignment-1))+SPIN_INPUT_OFFSET')
    source=replace_once(source,
        '(reinterpret_cast<std::uintptr_t>(guarded.data()+4)+63)&~std::uintptr_t(63)',
        '(reinterpret_cast<std::uintptr_t>(guarded.data()+4)+memoryAlignment-1)&~std::uintptr_t(memoryAlignment-1)')
    source=replace_once(source,
        'k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);',
        '''const auto inputAdvice=k::workspace_routing::adviseOwned(inputBase,n*16);
    const auto scratchAdvice=k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);''')
    source=replace_once(source,'<<" input_mod64="<<(reinterpret_cast<std::uintptr_t>(inputBase)%64)',
        '''<<" input_mod64="<<(reinterpret_cast<std::uintptr_t>(inputBase)%64)
             <<" owned_alignment="<<memoryAlignment<<" input_thp="<<inputAdvice.collapsed
             <<" scratch_thp="<<scratchAdvice.collapsed''')
    return source


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    args=parser.parse_args()
    print(generate(args.input))
