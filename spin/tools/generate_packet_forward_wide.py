"""Generate the K16 256-bit forward kernel from the retained exact schedule.

Keep the shared XOR circuit in 256-bit words. GFNI still operates on the two
packed 128-bit planes, and scratch retains its public internal planar layout.
This removes the per-lane output staging and shares route/address work without
enlarging the outer working set. See experiments/wide-k16-20261002 for receipts.
"""
import re

HELPERS = r'''
namespace wide256k16 {
struct V {__m512i lo,hi;};
static SPIN_FORCEINLINE __m256i vx(__m256i a,__m256i b){return _mm256_xor_si256(a,b);}
static SPIN_FORCEINLINE V px(V a,V b){return {_mm512_xor_si512(a.lo,b.lo),_mm512_xor_si512(a.hi,b.hi)};}
static SPIN_FORCEINLINE V pt(V a,V b,V c,int){return {_mm512_ternarylogic_epi64(a.lo,b.lo,c.lo,0x96),_mm512_ternarylogic_epi64(a.hi,b.hi,c.hi,0x96)};}
static SPIN_FORCEINLINE V join(__m256i a,__m256i b,__m256i c,__m256i d){return {_mm512_inserti64x4(_mm512_castsi256_si512(a),b,1),_mm512_inserti64x4(_mm512_castsi256_si512(c),d,1)};}
static SPIN_FORCEINLINE V broadcast(__m256i a){const auto x=_mm512_broadcast_i64x4(a);return {x,x};}
static SPIN_FORCEINLINE V load(const Block* raw,std::uint32_t offset,std::size_t stride){
 const auto a=_mm512_loadu_si512(raw+offset),b=_mm512_loadu_si512(raw+stride+offset);
 return {_mm512_permutex2var_epi64(a,_mm512_setr_epi64(0,1,8,9,2,3,10,11),b),_mm512_permutex2var_epi64(a,_mm512_setr_epi64(4,5,12,13,6,7,14,15),b)};
}
template<unsigned Mask,unsigned Degree> static SPIN_FORCEINLINE void storeLowMoments(V x,__m256i* moments){
 const auto folded=_mm512_xor_si512(x.lo,x.hi);const auto odd=_mm512_extracti64x4_epi64(folded,1);
 moments[Mask]=vx(_mm512_castsi512_si256(folded),odd);
 if constexpr(Degree>=1){const auto last=_mm512_extracti64x4_epi64(x.hi,1);moments[Mask+1]=odd;moments[Mask+2]=vx(_mm512_castsi512_si256(x.hi),last);if constexpr(Degree==2)moments[Mask+3]=last;}
}
template<bool Stream> struct Output {Block* out;SPIN_FORCEINLINE void operator()(std::size_t p,V v){
 if constexpr(Stream){_mm512_stream_si512(reinterpret_cast<__m512i*>(out+8*p),v.lo);_mm512_stream_si512(reinterpret_cast<__m512i*>(out+8*p+4),v.hi);}
 else {_mm512_storeu_si512(out+8*p,v.lo);_mm512_storeu_si512(out+8*p+4,v.hi);}
}};
'''

RUN = r'''
template<bool Stream> static SPIN_NOINLINE void run(const Block* input,Block* output,const Plan& p,const ForwardPlan& f){
 PackedState states[2]{},feedback;PackedExtra extras[2]{};
 alignas(64) __m128i unpacked[2][16],extraWords[2][4],syndrome128[16],extraSyndrome[4];
 alignas(64) __m256i words[16],ext[4],moments[64],syndrome[16];
 alignas(64) V packets[16];Output<Stream> emit{output};
 for(std::size_t epoch=0;epoch<p.n/64;++epoch){
  const auto* addresses=p.route.data()+16*epoch;const bool first=epoch==0;
  if(first){for(unsigned h=0;h<16;++h){packets[h]=load(input,addresses[h],p.scratchBlocks());emit(h,packets[h]);}}
  else {
   for(unsigned lane=0;lane<2;++lane){wideUnpackVbmi(states[lane],unpacked[lane]);unpackExtra<4>(extras[lane],extraWords[lane]);}
   for(unsigned j=0;j<16;++j)words[j]=_mm256_inserti128_si256(_mm256_castsi128_si256(unpacked[0][j]),unpacked[1][j],1);
   for(unsigned j=0;j<4;++j)ext[j]=_mm256_inserti128_si256(_mm256_castsi128_si256(extraWords[0][j]),extraWords[1][j],1);
   streamStepBorder<4>(words,ext,input,addresses,16*epoch,p.scratchBlocks(),emit,moments);
  }
  if(epoch+1==p.n/64)continue;
  if(first)packetMoments(packets,moments);
  finish(moments,syndrome);
  const auto finishLane=[&]<unsigned Lane>(){
   for(unsigned j=0;j<16;++j)syndrome128[j]=_mm256_extracti128_si256(syndrome[j],Lane);
   extraSyndrome[0]=_mm256_extracti128_si256(moments[3],Lane);
   extraSyndrome[1]=_mm256_extracti128_si256(moments[5],Lane);
   extraSyndrome[2]=_mm256_extracti128_si256(moments[9],Lane);
   extraSyndrome[3]=_mm256_extracti128_si256(moments[17],Lane);
   widePackVbmi(syndrome128,feedback);const auto ef=packExtra<4>(extraSyndrome);
   if(first){states[Lane]=feedback;extras[Lane]=ef;}else update(states[Lane],extras[Lane],feedback,ef,f.updates[epoch]);
  };
  finishLane.template operator()<0>();finishLane.template operator()<1>();
 }
 _mm_sfence();
}
}
'''

def outer_k16(source):
    """Use immediate shuffles for K16; retain the old loader at other sizes."""
    source = source.replace('template<unsigned L> static SPIN_FORCEINLINE __m512i loadQuad(const Block* in) {',
        """template<unsigned L,unsigned Lane,bool PairLoads> static SPIN_FORCEINLINE __m512i loadQuad(const Block* in) {
    if constexpr(L==2 && PairLoads)
        return _mm512_shuffle_i32x4(_mm512_loadu_si512(in),_mm512_loadu_si512(in+4),Lane==0?0x88:0xdd);
    in+=Lane;""")
    source = source.replace('template<unsigned L,unsigned S,unsigned P>', 'template<unsigned L,unsigned Lane,bool PairLoads,unsigned S,unsigned P>')
    source = source.replace('loadQuad<L>', 'loadQuad<L,Lane,PairLoads>')
    source = source.replace('template<unsigned L> static SPIN_NOINLINE void groupForward', 'template<unsigned L,unsigned Lane,bool PairLoads=false> static SPIN_NOINLINE void groupForward')
    source = source.replace('packMessage<L,', 'packMessage<L,Lane,PairLoads,')
    source = source.replace('groupForward<L>(in+256*g*L,', 'groupForward<L,0>(in+256*g*L,')
    for lane in (1,2,3):
        source = source.replace(f'groupForward<L>(in+256*g*L+{lane},', f'groupForward<L,{lane}>(in+256*g*L,')
    first = '        groupForward<L,0>(in+256*g*L,out+groupStride*g,f.outer.data()+16*g,nullptr);'
    source = source.replace(first, """        // Restrict the tuned loader to the measured K16 workload.
        if constexpr(L==2) {
            if(p.k==65536) {
                groupForward<2,0,true>(in+512*g,out+groupStride*g,f.outer.data()+16*g,nullptr);
                groupForward<2,1,true>(in+512*g,out+p.scratchBlocks()+groupStride*g,f.outer.data()+16*g,nullptr);
                continue;
            }
        }
""" + first)
    return source

def inner_k16(source, large):
    """Lift the unrolled schedule, without changing the underlying linear map."""
    def translate(text):
        return (text.replace('__m128i','__m256i').replace('__m512i','V')
            .replace('_mm512_xor_si512','px').replace('_mm512_ternarylogic_epi64','pt')
            .replace('_mm512_broadcast_i32x4','broadcast'))
    moments = large[large.index('struct FourMoments'):large.index('template<unsigned Mask')]
    moments += large[large.index('static SPIN_FORCEINLINE void packetMoments('):large.index('template<std::size_t... I>')]
    moments = translate(moments)
    moments = re.sub(r'_mm512_loadu_si512\(packets \+ (First(?: \+ \d)?)\)',r'packets[\1]',moments)
    stream = source[source.index('template<unsigned Extra,class Emit>'):source.index('\n}\n\nnamespace spin::detail::packet {')]
    stream = translate(stream).replace('const block* raw','const Block* raw')
    stream = stream.replace('std::size_t packetBase,Emit& emit','std::size_t packetBase,std::size_t stride,Emit& emit')
    stream = re.sub(r'_mm512_loadu_si512\(raw\+route\[(\d+)\]\)',r'load(raw,route[\1],stride)',stream)
    stream = stream.replace('route,packetBase,emit,high','route,packetBase,stride,emit,high')
    where = source.index('void forwardInnerWide(')
    source = source[:where]+HELPERS+moments+stream+RUN+source[where:]
    source = source.replace('if(lanes==2)innerWide<2>(in,out,p,f);else innerWide<4>(in,out,p,f);',
        """if(lanes==2 && p.k==65536) {
        if((reinterpret_cast<std::uintptr_t>(out)&63)==0)wide256k16::run<true>(in,out,p,f);
        else wide256k16::run<false>(in,out,p,f);
    } else if(lanes==2)innerWide<2>(in,out,p,f);else innerWide<4>(in,out,p,f);""")
    return source
