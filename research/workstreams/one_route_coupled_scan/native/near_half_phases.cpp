// Measurement-only driver. Reuse the frozen harness and its setup verbatim;
// rename its entry point so this translation unit has one diagnostic main.
#define SPIN_OR_WINDOW_NINE 1
#define SPIN_OR_WINDOW_COORDINATES 1
#define SPIN_OR_WINDOW_RAW_SCALE 1
#define main retained_near_half_main
#include "near_half_native.cpp"
#undef main
#include "NearHalfPhaseKernel.h"

namespace phases {
using namespace native;
namespace pk=spin::research::one_route::phase;
struct Times {
    double outer=0,spectra=0,inverse=0,route=0,append=0,total=0;
};
double milliseconds(Clock::time_point a,Clock::time_point b){
    return std::chrono::duration<double,std::milli>(b-a).count();
}
NI void productPhase(const Kernel::State& s,const Kernel::PreparedBatch& k,Kernel::Scratch& w){
    pk::spectra(s,k,w);
}
NI void inversePhase(Kernel::Scratch& w,E* feedback){pk::finish(w,feedback);}
NI void appendPhase(Kernel::State& s,const E* output){s.append(output);}
NI void routePhase(const Plan& p,const E* input,E* output,const E* feedback,
                   std::size_t offset,std::size_t base){
    if(p.blocks>=256){
        const auto* source=input+base;const auto* labels=p.route.data()+offset;
        for(unsigned j=0;j<256;j+=4){
            auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j])));
            x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+1])),1);
            x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+2])),2);
            x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[j+3])),3);
            _mm512_store_si512(output+offset+j,_mm512_xor_si512(x,_mm512_load_si512(feedback+j)));
        }
    }else{
        for(unsigned j=0;j<256;++j){
            const auto index=offset+j,region=(index/p.blocks)*p.blocks;
            output[index]=input[region+p.route[index]]^feedback[j];
        }
    }
}
template<bool Timed> NI Times splitFull(const Plan& p,const E* input,E* workspace,E* output){
    Times result;
    const auto beginning=Clock::now();
    outer(p,input,workspace);
    auto mark=Clock::now();
    if constexpr(Timed)result.outer=milliseconds(beginning,mark);
    Kernel::State state;state.reset();Kernel::Scratch scratch;
    alignas(64) E feedback[256];
    std::size_t base=0,next=p.blocks;
    for(std::size_t t=0;t<p.batches;++t){
        if constexpr(Timed)mark=Clock::now();
        productPhase(state,p.keys[t],scratch);
        if constexpr(Timed){auto now=Clock::now();result.spectra+=milliseconds(mark,now);mark=now;}
        inversePhase(scratch,feedback);
        if constexpr(Timed){auto now=Clock::now();result.inverse+=milliseconds(mark,now);mark=now;}
        const auto offset=256*t;
        if(offset==next){base=next;next+=p.blocks;}
        routePhase(p,workspace,output,feedback,offset,base);
        if constexpr(Timed){auto now=Clock::now();result.route+=milliseconds(mark,now);mark=now;}
        appendPhase(state,output+offset);
        if constexpr(Timed){auto now=Clock::now();result.append+=milliseconds(mark,now);}
    }
    result.total=milliseconds(beginning,Clock::now());
    return result;
}
template<bool HotKeys> NI void spectraStream(const Kernel::State& s,const Plan& p,Kernel::Scratch& w){
    for(std::size_t t=0;t<p.batches;++t){
        pk::spectra(s,p.keys[HotKeys?0:t],w);
        // Observe every iteration, including repeated equal hot-key results.
        asm volatile("" : : "g"(w.spectra) : "memory");
    }
}
// Output layout diagnostic, not a replacement encoder. Both modes evaluate
// the identical four-block RS circuit and read identical input. Only stores
// differ: coordinate-major plane layout versus contiguous four-block tiles.
template<bool Planar,bool Arithmetic> NI void outerDiagnostic(const Plan& p,const E* input,E* output){
    alignas(64) E scratch[256*4];
    for(std::size_t row=0;row<p.blocks;row+=4){
        for(unsigned j=0;j<128;++j){
            scratch[4*j]=input[(row+0)*128+j];
            scratch[4*j+1]=input[(row+1)*128+j];
            scratch[4*j+2]=input[(row+2)*128+j];
            scratch[4*j+3]=input[(row+3)*128+j];
        }
        std::memset(scratch+512,0,512*sizeof(E));
        if constexpr(Arithmetic)FFT::evaluateNovelFour(scratch);
        for(unsigned j=0;j<256;++j){
            auto* dest=Planar?output+j*p.blocks+row:output+256*row+4*j;
            _mm512_store_si512(dest,_mm512_load_si512(scratch+4*j));
        }
    }
}
double median(std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];}
void run(unsigned logK,unsigned repeats){
    require(logK>=12&&logK<=20&&repeats>=3&&repeats<=101,"bounded phase args");
    Plan p(std::size_t(1)<<logK);
    Aligned<E> input(p.K),workspace(p.N),output(p.N),reference(p.N);
    std::mt19937_64 rng(931);for(std::size_t i=0;i<p.K;++i)input[i]=random(rng);
    full<false,true>(p,input.data(),workspace.data(),reference.data());
    const auto expected=checksum(reference.data(),p.N);
    splitFull<false>(p,input.data(),workspace.data(),output.data());
    require(std::memcmp(output.data(),reference.data(),p.N*sizeof(E))==0,"untimed split changed encoder output");
    splitFull<true>(p,input.data(),workspace.data(),output.data());
    require(std::memcmp(output.data(),reference.data(),p.N*sizeof(E))==0,"timed split changed encoder output");
    std::cout<<std::fixed<<std::setprecision(6)<<"K="<<p.K<<" N="<<p.N<<" repeats="<<repeats
        <<" checksum="<<expected<<" split_oracle=PASS\n"<<std::flush;
    std::vector<double> original,untimed,total,ot,sp,iv,rt,ap,other;
    for(unsigned trial=0;trial<repeats+3;++trial){
        // Alternate order to avoid always assigning the same temperature/order
        // effects to one mode. No experiments execute concurrently.
        double a,b;Times c;
        auto retained=[&]{auto start=Clock::now();full<false,true>(p,input.data(),workspace.data(),output.data());return milliseconds(start,Clock::now());};
        if(trial&1){c=splitFull<true>(p,input.data(),workspace.data(),output.data());b=splitFull<false>(p,input.data(),workspace.data(),output.data()).total;a=retained();}
        else{a=retained();b=splitFull<false>(p,input.data(),workspace.data(),output.data()).total;c=splitFull<true>(p,input.data(),workspace.data(),output.data());}
        require(checksum(output.data(),p.N)==expected,"timing output changed");
        if(trial<3)continue;
        original.push_back(a);untimed.push_back(b);total.push_back(c.total);ot.push_back(c.outer);
        sp.push_back(c.spectra);iv.push_back(c.inverse);rt.push_back(c.route);ap.push_back(c.append);
        other.push_back(c.total-c.outer-c.spectra-c.inverse-c.route-c.append);
    }
    const auto show=[](const char* name,double ms){std::cout<<name<<"_ms="<<ms<<"\n"<<std::flush;};
    show("retained_full",median(original));show("split_untimed",median(untimed));show("split_timed",median(total));
    show("outer",median(ot));show("mix_products_reduction_keys",median(sp));show("inverse_transforms_extract",median(iv));
    show("routed_reads_xor_output",median(rt));show("state_append_forward_transforms",median(ap));show("unattributed_timer_loop",median(other));

    // Key reuse changes the setup and is ONLY a cache diagnostic. This does
    // not produce encoded words and is not a candidate certified construction.
    Kernel::State state;state.load(workspace.data());Kernel::Scratch scratch;
    show("spectra_streamed_keys",measure(repeats,[&]{spectraStream<false>(state,p,scratch);}));
    std::cout<<"streamed_checksum="<<(checksum(scratch.spectra[0],256)^checksum(scratch.spectra[1],256))<<"\n";
    show("spectra_one_hot_key",measure(repeats,[&]{spectraStream<true>(state,p,scratch);}));
    std::cout<<"hot_checksum="<<(checksum(scratch.spectra[0],256)^checksum(scratch.spectra[1],256))<<"\n";
    show("prepared_key_read",measure(repeats,[&]{readTape(p.keys.data(),p.batches,output.data());}));

    show("outer_plane_layout",measure(repeats,[&]{outerDiagnostic<true,true>(p,input.data(),output.data());}));
    outer(p,input.data(),workspace.data());
    require(std::memcmp(output.data(),workspace.data(),p.N*sizeof(E))==0,"outer diagnostic planar mismatch");
    show("outer_tiled_layout",measure(repeats,[&]{outerDiagnostic<false,true>(p,input.data(),output.data());}));
    for(std::size_t row=0;row<p.blocks;row+=4)for(unsigned j=0;j<256;++j)
        require(std::memcmp(output.data()+row*256+4*j,workspace.data()+j*p.blocks+row,4*sizeof(E))==0,"outer diagnostic tiled mismatch");
    show("outer_plane_no_fft",measure(repeats,[&]{outerDiagnostic<true,false>(p,input.data(),output.data());}));
    show("outer_tiled_no_fft",measure(repeats,[&]{outerDiagnostic<false,false>(p,input.data(),output.data());}));
    std::cout<<"layout_checks=PASS final_checksum="<<checksum(output.data(),p.N)<<"\n";
}
} //namespace phases
int main(int argc,char** argv){
    try{
        if(argc!=3){std::cerr<<"usage: near_half_phases logK repeats\n";return 2;}
        phases::run(unsigned(std::stoul(argv[1])),unsigned(std::stoul(argv[2])));return 0;
    }catch(const std::exception& e){std::cerr<<"FAIL: "<<e.what()<<"\n";return 1;}
}
