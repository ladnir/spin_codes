// Exact deterministic reformulation; sampled geometry is NOT an ensemble proof.
// Serial, no files, <256MiB owned arrays; remote runner limits time and memory.
#define main quartet_original_diagnostic_main
#include "quartet_support_cut.cpp"
#undef main
#include <array>

struct Window {
    std::uint32_t first;
    std::uint8_t length,weight,degree;
    std::array<std::uint8_t,17> assigned{};
};
static_assert(sizeof(Window)==24);

struct Allocation {
    const Geometry &g;
    std::vector<std::uint16_t> quartet;
    std::vector<Window> windows;
    std::vector<int> load;
    std::vector<std::uint32_t> order;
    i64 empty=0,total=0;
    explicit Allocation(const Geometry &geometry):g(geometry),quartet(g.N),load(g.N,0) {
        windows.reserve(2*g.N+3*g.width);
        for(int t=0;t<g.steps;++t)for(int base=0;base<g.width;base+=4)
            for(int k=0;k<4;++k)quartet[t*g.width+g.groups[t*g.width+base+k]]=base;
        auto append=[&](int t,int j,int length,int weight) {
            Window w{};w.first=t*g.width+j;w.length=length;w.weight=weight;
            int labels[17];w.degree=neighbors(w,labels);windows.push_back(w);
            total+=weight;if(!w.degree)empty+=weight;
        };
        assert(g.steps>=4);
        for(int t=0;t<g.steps;++t)for(int j=0;j<g.width;++j) {
            if(t+3<=g.steps)append(t,j,3,1);
            if(t+4<=g.steps)append(t,j,4,24);
        }
        for(int j=0;j<g.width;++j) {
            append(g.steps-1,j,1,25);append(g.steps-2,j,2,25);append(g.steps-3,j,3,24);
        }
        assert(total==25LL*g.N);
        std::array<std::size_t,19> begin{};
        for(const auto &w:windows)++begin[w.degree+1];
        for(int i=1;i<19;++i)begin[i]+=begin[i-1];
        auto next=begin;order.resize(windows.size());
        for(std::uint32_t i=0;i<windows.size();++i)order[next[windows[i].degree]++]=i;
        std::size_t bytes=quartet.capacity()*2+windows.capacity()*sizeof(Window)
            +load.capacity()*4+order.capacity()*4+std::size_t(g.N)*12;
        assert(bytes<256ULL*(1<<20));
    }
    int neighbors(const Window &w,int *out) const {
        int count=0,t=w.first/g.width,j=w.first%g.width;
        if(g.route[w.first]<g.inactive)out[count++]=w.first;
        for(int k=1;k<=w.length && t+k<g.steps;++k) {
            int offset=(t+k)*g.width,base=quartet[offset+j];
            for(int b=0;b<4;++b) {
                int p=offset+g.groups[offset+base+b];
                if(g.route[p]<g.inactive)out[count++]=p;
            }
        }
        assert(count<=17);return count;
    }
    // Every mass unit is assigned to an actual blocking label. Trimming
    // overloaded labels gives a feasible dual certificate, regardless of order.
    void pass(bool reverse) {
        for(std::size_t pos=0;pos<order.size();++pos) {
            auto &w=windows[order[reverse?order.size()-1-pos:pos]];
            if(!w.degree)continue;
            int labels[17];assert(neighbors(w,labels)==w.degree);
            for(int k=0;k<w.degree;++k){load[labels[k]]-=w.assigned[k];assert(load[labels[k]]>=0);w.assigned[k]=0;}
            for(int unit=0;unit<w.weight;++unit) {
                int best=0;
                for(int k=1;k<w.degree;++k)if(load[labels[k]]<load[labels[best]])best=k;
                ++load[labels[best]];++w.assigned[best];
            }
        }
    }
    i64 bound(bool verify=true) const {
        i64 result=empty,mass=empty;
        for(int value:load){assert(value>=0);result+=std::max(0,value-100);mass+=value;}
        assert(mass==total);
        if(verify) {
            std::vector<int> replay(g.N,0);
            for(const auto &w:windows) {
                int labels[17],count=neighbors(w,labels),sum=0;
                assert(count==w.degree);
                for(int k=0;k<count;++k){replay[labels[k]]+=w.assigned[k];sum+=w.assigned[k];}
                assert(sum==(count?w.weight:0));
            }
            assert(replay==load);
        }
        return result;
    }
    void uniform_bound() const {
        constexpr i64 scale=12252240;
        std::vector<i64> uniform(g.N,0);
        for(const auto &w:windows)if(w.degree) {
            int labels[17];assert(neighbors(w,labels)==w.degree);assert(scale%w.degree==0);
            i64 share=scale*w.weight/w.degree;
            for(int k=0;k<w.degree;++k)uniform[labels[k]]+=share;
        }
        i64 excess=0;
        for(i64 value:uniform)excess+=std::max(i64(0),value-100*scale);
        std::cout<<"uniform_bound_numerator="<<empty*scale+excess<<" denominator="<<100*scale
                 <<" uniform_bound_float="<<double(empty*scale+excess)/(100*scale)<<std::endl;
    }
};

void window_tests() {
    i64 histories=0,labelsets=0;
    for(int seed=0;seed<8;++seed) {
        Geometry g(4,8,2,2,seed);Allocation a(g);
        std::vector<int> inactive;
        for(int p=0;p<g.N;++p)if(g.route[p]<g.inactive)inactive.push_back(p);
        std::vector<i64> best(1u<<inactive.size(),std::numeric_limits<i64>::min()/4);
        i64 optimum=0;
        std::vector<unsigned char> state(g.N);
        for(unsigned mask=0;mask<(1u<<g.N);++mask) {
            for(int p=0;p<g.N;++p)state[p]=(mask>>p)&1;
            auto c=count(g,state);i64 score=25*c.A-100*c.C-74*c.D;optimum=std::max(optimum,score);
            unsigned covered=0;
            for(unsigned i=0;i<inactive.size();++i) {
                int p=inactive[i],t=p/g.width,j=p%g.width;
                bool hit=state[p];
                if(t)for(int k=0;k<4;++k)hit|=state[(t-1)*g.width+g.groups[t*g.width+a.quartet[p]+k]];
                if(hit)covered|=1u<<i;
            }
            assert(std::popcount(covered)==c.C);
            best[covered]=std::max(best[covered],25*c.A-74*c.D);++histories;
        }
        i64 window_optimum=0;
        for(unsigned paid=0;paid<best.size();++paid) {
            unsigned cleared=0;
            for(unsigned i=0;i<inactive.size();++i)if((paid>>i)&1)cleared|=1u<<inactive[i];
            i64 expected=0;
            for(unsigned subset=paid;;subset=(subset-1)&paid){expected=std::max(expected,best[subset]);if(!subset)break;}
            i64 value=0;
            for(const auto &w:a.windows) {
                int labels[17],n=a.neighbors(w,labels);bool enabled=true;
                for(int k=0;k<n;++k)enabled&=(cleared>>labels[k])&1;
                if(enabled)value+=w.weight;
            }
            assert(value==expected);
            window_optimum=std::max(window_optimum,value-100*std::popcount(paid));++labelsets;
        }
        assert(window_optimum==optimum);
        a.pass(false);assert(optimum<=a.bound());a.pass(true);assert(optimum<=a.bound());
    }
    std::cout<<"exhaustive_window_history_checks="<<histories<<" exact_paid_label_sets="<<labelsets
             <<" free_terminal_identity=1 assignment_certificates_verified=1"<<std::endl;
}

int main(int argc,char **argv) {
    auto start=std::chrono::steady_clock::now();window_tests();
    if(argc>1 && std::string(argv[1])=="--tests-only")return 0;
    Geometry g(512,8192,256,2100,20261027);Allocation a(g);
    std::cout<<"empty_window_weight="<<a.empty<<" denominator=100 total_weight="<<a.total
             <<" window_count="<<a.windows.size()<<" exact_for_all_supports_on_this_geometry=1 ensemble_tail=0"<<std::endl;
    a.uniform_bound();
    for(int pass=0;pass<4;++pass) {
        a.pass(pass%2);auto upper=a.bound();
        std::cout<<"allocation_pass="<<pass+1<<" bound_numerator="<<upper<<" denominator=100 bound_float="<<double(upper)/100
                 <<" known_sample_optimum=3988.44 sufficient_ensemble_offset=18171 ensemble_tail=0"<<std::endl;
    }
    std::cout<<"seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<" files_written=0"<<std::endl;
}
