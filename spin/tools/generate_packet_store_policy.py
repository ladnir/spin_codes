"""Apply per-call output store selection to the retained forward schedules.

Keep this final transformation after deriving wide and gathered circuits so
their data flow stays independent of the public output-store policy.
"""
import re

STORE = '''template<bool Stream> static SPIN_FORCEINLINE void storeOutput(Block* out,__m512i value) {
    if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(out),value);
    else _mm512_storeu_si512(out,value);
}
'''

def replace_once(source, before, after):
    assert source.count(before) == 1, 'store-policy generator anchor changed: ' + before[:80]
    return source.replace(before, after)

def outer(source):
    source = replace_once(source, 'static SPIN_FORCEINLINE void unpackHalf(',
                          STORE + 'template<bool Stream> static SPIN_FORCEINLINE void unpackHalf(')
    source = re.sub(r'_mm512_stream_si512\(reinterpret_cast<__m512i\*>\((out\+route\[\d+\])\),',
                    r'storeOutput<Stream>(\1,', source)
    # The newly inserted store helper is not affected by this pointer-specific
    # rewrite. It is the only NT instruction in the final outer source.
    assert source.count('_mm512_stream_si512') == 1
    source = replace_once(source, 'static SPIN_NOINLINE void groupForward(',
                          'template<bool Stream> static SPIN_NOINLINE void groupForward(')
    source = source.replace('unpackHalf(v,', 'unpackHalf<Stream>(v,')
    source = source.replace('unpackHalf(v+4,', 'unpackHalf<Stream>(v+4,')
    start = source.index('void forwardOuter(')
    source = source[:start] + '''void forwardOuter(const Block* in,Block* out,const Plan& p,const ForwardPlan& f,bool stream) {
    // For large K this can be caller output reused as the inner's input.
    // Cached applies to those intermediate writes as well as final writes.
    if(stream && (reinterpret_cast<std::uintptr_t>(out)&63U)==0) {
        for(std::size_t g=0;g<p.groups;++g)
            groupForward<true>(in+256*g,out,f.outer.data()+16*g,f.inverseRoute.data()+groupStride*g/4);
        _mm_sfence();
    } else {
        for(std::size_t g=0;g<p.groups;++g)
            groupForward<false>(in+256*g,out,f.outer.data()+16*g,f.inverseRoute.data()+groupStride*g/4);
    }
}
}
'''
    return source

def inner(source, name):
    start = source.index('struct DirectOutput {')
    end = source.index('template<unsigned Extra,class Emit>', start)
    source = source[:start] + STORE + '''template<bool Stream> struct DirectOutput {
    Block* output;
    SPIN_FORCEINLINE void operator()(std::size_t p,__m512i value) {storeOutput<Stream>(output+4*p,value);}
};
''' + source[end:]
    start = source.index('void '+name+'(')
    source = source[:start] + f'''void {name}(const Block* in,Block* out,const Plan& p,const ForwardPlan& f,bool stream) {{
    if(stream && (reinterpret_cast<std::uintptr_t>(out)&63U)==0) {{
        DirectOutput<true> emit{{out}};
        forwardBorder<4>(in,p.n,f.updates.data(),p.route.data(),emit);
        _mm_sfence();
    }} else {{
        DirectOutput<false> emit{{out}};
        forwardBorder<4>(in,p.n,f.updates.data(),p.route.data(),emit);
    }}
}}
}}
'''
    return source

def wide(source):
    source = replace_once(source,
        'template<unsigned L> void innerWide(const Block* in,Block* out,const Plan& p,const ForwardPlan& f) {\n'
        '    if((reinterpret_cast<std::uintptr_t>(out)&63)==0)\n'
        '        forwardBorder<L,true>(in,p.n,f.updates.data(),p.route.data(),out,p.scratchBlocks());\n'
        '    else forwardBorder<L,false>(in,p.n,f.updates.data(),p.route.data(),out,p.scratchBlocks());\n'
        '    _mm_sfence();\n}',
        '''template<unsigned L> void innerWide(const Block* in,Block* out,const Plan& p,const ForwardPlan& f,bool stream) {
    if(stream && (reinterpret_cast<std::uintptr_t>(out)&63U)==0) {
        forwardBorder<L,true>(in,p.n,f.updates.data(),p.route.data(),out,p.scratchBlocks());
        _mm_sfence();
    } else forwardBorder<L,false>(in,p.n,f.updates.data(),p.route.data(),out,p.scratchBlocks());
}''')
    # The tuned K16 kernel is already store-specialized; make its fence equally
    # specialized and pass the public policy into its once-per-call dispatcher.
    source = replace_once(source, '\n _mm_sfence();\n', '\n if constexpr(Stream)_mm_sfence();\n')
    source = replace_once(source,
        'void forwardInnerWide(const Block* in,Block* out,const Plan& p,const ForwardPlan& f,unsigned lanes) {',
        'void forwardInnerWide(const Block* in,Block* out,const Plan& p,const ForwardPlan& f,unsigned lanes,bool stream) {')
    source = replace_once(source,
        'if((reinterpret_cast<std::uintptr_t>(out)&63)==0)wide256k16::run<true>',
        'if(stream && (reinterpret_cast<std::uintptr_t>(out)&63U)==0)wide256k16::run<true>')
    source = source.replace('innerWide<2>(in,out,p,f)', 'innerWide<2>(in,out,p,f,stream)')
    source = source.replace('innerWide<4>(in,out,p,f)', 'innerWide<4>(in,out,p,f,stream)')
    return source

def apply_output_stores(name, source):
    if name == 'Outer': return outer(source)
    if name == 'Inner': return inner(source, 'forwardInner')
    if name == 'InnerGather': return inner(source, 'forwardInnerGather')
    if name == 'InnerWide': return wide(source)
    return source
