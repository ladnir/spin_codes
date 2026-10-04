#pragma once
#include <spin/Code.h>
#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <span>
#include <vector>

namespace spin::test {
template<unsigned Bytes> struct alignas(16) StoreRecord {std::array<std::byte,Bytes> bytes;};
template<unsigned Bytes>
void typedForward(const Code& code,std::span<const std::byte> input,std::span<std::byte> output,
                  Workspace& work,OutputStores stores) {
    using Record=StoreRecord<Bytes>;
    code.forward<Record>({reinterpret_cast<const Record*>(input.data()),input.size()/Bytes},
        {reinterpret_cast<Record*>(output.data()),output.size()/Bytes},work,stores);
}
inline void typedForward(const Code& code,std::span<const std::byte> input,std::span<std::byte> output,
                         Workspace& work,OutputStores stores) {
    switch(work.width()) {
    case Width::Bits128:typedForward<16>(code,input,output,work,stores);break;
    case Width::Bits256:typedForward<32>(code,input,output,work,stores);break;
    case Width::Bits512:typedForward<64>(code,input,output,work,stores);break;
    }
}

// Exercise the public call boundary, not private kernel selectors. Every policy
// must return completed output that is immediately readable by the caller.
// A missing streaming-store fence can therefore fail the immediate comparison.
template<class AllocationCheck>
void forwardOutputStores(const Code& code,AllocationCheck allocationCheck) {
    const auto descriptor=code.descriptor();
    for(auto width:{Width::Bits128,Width::Bits256,Width::Bits512}) {
        if(!code.supports_forward(width))continue;
        const std::size_t bytes=static_cast<unsigned>(width),inBytes=code.message_size()*bytes,outBytes=code.code_size()*bytes;
        auto work=code.make_workspace(width);const auto scratchBytes=work.bytes();
        Buffer input(inBytes+64),output(outBytes+128),expected(outBytes);
        auto in=input.bytes().subspan(16,inBytes);
        std::vector<std::byte> before(inBytes);
        unsigned pattern=0,offset=0;OutputStores selected=OutputStores::Streaming;
        const auto check=[&](bool pass,const char* message) {
            if(pass)return;
            std::cerr<<"output-store test: family="<<static_cast<unsigned>(code.specification().parameters)
                <<" K="<<code.message_size()<<" width_bits="<<bytes*8
                <<" backend="<<static_cast<unsigned>(code.backend())<<" pattern="<<pattern
                <<" alignment_offset="<<offset<<" policy="<<static_cast<unsigned>(selected)<<'\n';
            throw std::runtime_error(message);
        };
        const auto matches=[&](std::span<const std::byte> actual,std::span<const std::byte> reference) {
            check(actual.size()==reference.size(),"store-policy comparison shape");
            if(!std::memcmp(actual.data(),reference.data(),actual.size()))return;
            std::size_t first=0;while(actual[first]==reference[first])++first;
            std::cerr<<"first different byte="<<first<<" record="<<first/bytes<<'\n';
            check(false,"forward store policy changed output");
        };
        const auto verify=[&](std::span<const std::byte> out) {
            matches(out,expected.bytes());matches(in,before);
            check(std::all_of(output.bytes().begin(),output.bytes().begin()+64+offset,[](auto x){return x==std::byte{0xa5};}),"store-policy prefix guard");
            check(std::all_of(output.bytes().begin()+64+offset+outBytes,output.bytes().end(),[](auto x){return x==std::byte{0xa5};}),"store-policy suffix guard");
            check(code.descriptor()==descriptor && work.bytes()==scratchBytes,"store policy changed descriptor or workspace");
        };
        for(pattern=0;pattern<3;++pattern) {
            std::uint64_t seed=919+pattern;
            for(std::size_t i=0;i<in.size();i+=8) {
                seed^=seed<<13;seed^=seed>>7;seed^=seed<<17;
                std::memcpy(in.data()+i,&seed,8);
            }
            if(pattern)std::fill(in.begin(),in.end(),std::byte{});
            if(pattern==2) {in.front()=std::byte{1};in.back()=std::byte{128};}
            std::copy(in.begin(),in.end(),before.begin());
            allocationCheck([&]{code.forward_bytes(in,expected.bytes(),work);});
            // The same dirty workspace and output change policy repeatedly.
            // Streaming may fall back at 16-byte but not 64-byte alignment.
            for(auto alignment:{0U,16U}) {
                offset=alignment;auto out=output.bytes().subspan(64+offset,outBytes);
                for(auto policy:{OutputStores::Cached,OutputStores::Streaming,OutputStores::Cached}) {
                    selected=policy;std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
                    allocationCheck([&]{code.forward_bytes(in,out,work,selected);});
                    verify(out);
                }
                // Typed overrides must reach the same width-specific policy path.
                if(pattern==0)for(auto policy:{OutputStores::Cached,OutputStores::Streaming}) {
                    selected=policy;std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
                    allocationCheck([&]{typedForward(code,in,out,work,selected);});
                    verify(out);
                }
                selected=static_cast<OutputStores>(99);
                for(bool typed:{false,true}) {
                    std::fill(output.bytes().begin(),output.bytes().end(),std::byte{0xa5});
                    bool rejected=false;
                    try {
                        if(typed)typedForward(code,in,out,work,selected);
                        else code.forward_bytes(in,out,work,selected);
                    } catch(const std::invalid_argument&) {rejected=true;}
                    check(rejected,"unknown output-store policy accepted");
                    check(std::all_of(output.bytes().begin(),output.bytes().end(),[](auto x){return x==std::byte{0xa5};}),"invalid store policy wrote output");
                    matches(in,before);
                }
            }
        }
        // Changing the forward store preference cannot change transposed encoding.
        if(width==Width::Bits128) {
            selected=OutputStores::Cached;
            Buffer a(inBytes),b(inBytes);code.transpose_bytes(expected.bytes(),a.bytes(),work);
            allocationCheck([&]{code.forward_bytes(in,output.bytes().first(outBytes),work,OutputStores::Cached);});
            allocationCheck([&]{code.transpose_bytes(expected.bytes(),b.bytes(),work);});
            matches(a.bytes(),b.bytes());
        }
    }
}
}
