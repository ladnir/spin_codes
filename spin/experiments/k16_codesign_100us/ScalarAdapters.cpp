// Expose the retained independent literal outer, without editing its source.
#include "../../../research/workstreams/k16_design/implementation/RsScalar.cpp"
namespace spin::research::k16codesign {
void outerScalar(const rs::Block* input,rs::Block* output,const rs::Plan& plan) {
    for(std::size_t g=0;g<plan.groups;++g)
        rs::outerTranspose16(input+rs::groupStride*g,output+128*g,plan.outerMatrices16.data()+16*g);
}
void outerForwardScalar(const rs::Block* input,rs::Block* routed,const rs::Plan& plan) {
    for(std::size_t g=0;g<plan.groups;++g)
        rs::outerForward16(input+128*g,routed+rs::groupStride*g,plan.outerMatrices16.data()+16*g);
}
}
