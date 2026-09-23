    template<class E,class Op> static inline void imt(E* state,u32 u,u32 v,const Op& op) {
        E dot=op(state[0],state[0]);
        while(u) {const auto j=std::countr_zero(u);u&=u-1;dot=op(dot,state[j]);}
        while(v) {const auto j=std::countr_zero(v);v&=v-1;state[j]=op(state[j],dot);}
    }
    template<unsigned T,unsigned D,unsigned Base=0,class E,class Op> static inline void zetaStage(E* x,const Op& op) {
        // Same compile-time pruning and paired loop shape as the SIMD kernel.
        if constexpr(std::popcount(Base/(2*D))<=2) {
            if constexpr(D==1) x[Base]=op(x[Base],x[Base+1]);
            else for(unsigned j=0;j<D;j+=2) {
                x[Base+j]=op(x[Base+j],x[Base+j+D]);
                x[Base+j+1]=op(x[Base+j+1],x[Base+j+D+1]);
            }
        }
        if constexpr(Base+2*D<T) zetaStage<T,D,Base+2*D>(x,op);
    }
    template<unsigned T,unsigned D=T/2,class E,class Op> static inline void zeta(E* x,const Op& op) {
        zetaStage<T,D>(x,op);
        if constexpr(D>1) zeta<T,D/2>(x,op);
    }
    template<class Map,unsigned Rounds,ValueElement E,class Op,class Route> void run(const E* in,E* out,Workspace<E>& w,const Op& op,const Route& route,const u32* maskData) const {
        std::array<E,Map::S> state,syndrome;
        std::array<E,Map::T> raw;
        const auto epochs=code_size()/Map::T;
        auto emit=[&](std::size_t i,const E& v) {w.values_[route(i)]=v;};
        for(std::size_t epoch=epochs;epoch-->0;) {
            const auto base=epoch*Map::T;
            if constexpr(requires{route.begin_epoch(base);})route.begin_epoch(base);
            if(epoch+1==epochs) {
                for(unsigned j=Map::T;j-->0;) {raw[j]=in[base+j];emit(base+j,raw[j]);}
            } else Map::emitShared(in+base,raw.data(),state.data(),base,emit,op);
            if(epoch==0) break;
            zeta<Map::T>(raw.data(),op);Map::finish(raw.data(),syndrome.data(),op);
            if(epoch+1==epochs) state=syndrome;
            else {
                const auto* masks=maskData+2*Rounds*epoch;
                if constexpr(Rounds==2) imt(state.data(),masks[2],masks[3],op);
                imt(state.data(),masks[0],masks[1],op);
                for(unsigned j=0;j<Map::S;++j) state[j]=op(state[j],syndrome[j]);
            }
        }
        // All input was consumed above: exactly in-place output is now safe.
        for(std::size_t row=0;row<data_->k/128;++row)
            detail::generic::bch(w.values_.data()+256*row,out+128*row,op);
    }
