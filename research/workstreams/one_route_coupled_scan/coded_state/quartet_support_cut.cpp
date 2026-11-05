// Structural diagnostic, NOT an ensemble-distance certificate.
// One thread, memory bounded by explicit graph-size assertions. No files.
// Minimize C_inactive + (carry/denominator)*deaths - (numerator/denominator)*activity.
#include <algorithm>
#include <cassert>
#include <bit>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using i64 = std::int64_t;
struct Edge { int to, next, cap; };
struct Flow {
    std::vector<int> head, level, current, queue, initial;
    std::vector<Edge> edges;
    explicit Flow(int n, std::size_t reserve) : head(n,-1),level(n),current(n),queue(n) {
        edges.reserve(reserve); initial.reserve(reserve/2);
    }
    void add(int a,int b,int cap) {
        assert(cap>=0 && edges.size()+2<=edges.capacity());
        initial.push_back(cap);
        edges.push_back({b,head[a],cap}); head[a]=int(edges.size())-1;
        edges.push_back({a,head[b],0}); head[b]=int(edges.size())-1;
    }
    bool bfs(int source,int sink) {
        std::fill(level.begin(),level.end(),-1);int lo=0,hi=0;
        queue[hi++]=source;level[source]=0;
        while(lo<hi) {
            int v=queue[lo++];
            for(int e=head[v];e!=-1;e=edges[e].next) if(edges[e].cap && level[edges[e].to]<0) {
                level[edges[e].to]=level[v]+1;queue[hi++]=edges[e].to;
            }
        }
        return level[sink]>=0;
    }
    int dfs(int v,int sink,int amount) {
        if(v==sink)return amount;
        for(int &e=current[v];e!=-1;e=edges[e].next) {
            Edge &a=edges[e];
            if(a.cap && level[a.to]==level[v]+1) {
                int sent=dfs(a.to,sink,std::min(amount,a.cap));
                if(sent){a.cap-=sent;edges[e^1].cap+=sent;return sent;}
            }
        }
        return 0;
    }
    i64 solve(int source,int sink,bool verbose) {
        i64 value=0;int phases=0;
        while(bfs(source,sink)) {
            current=head;
            while(int sent=dfs(source,sink,std::numeric_limits<int>::max()))value+=sent;
            ++phases;
            if(verbose && (phases<=3 || phases%5==0))
                std::cout<<"flow_phase="<<phases<<" depth="<<level[sink]<<" value="<<value<<std::endl;
        }
        return value;
    }
    void verify(int source,int sink,i64 value) const {
        std::vector<i64> balance(head.size(),0);
        for(std::size_t e=0;e<edges.size();e+=2) {
            const auto &a=edges[e], &b=edges[e+1];
            assert(a.cap>=0 && b.cap>=0 && i64(a.cap)+b.cap==initial[e/2]);
            balance[b.to]-=b.cap;balance[a.to]+=b.cap;
        }
        for(std::size_t i=0;i<balance.size();++i)
            assert(balance[i]==(int(i)==source?-value:int(i)==sink?value:0));
    }
};

struct Geometry {
    int width,rows,regions,inactive,steps,N;
    std::vector<std::uint16_t> route,groups;
    Geometry(int w,int m,int b,int r,std::uint64_t seed):
        width(w),rows(m),regions(b),inactive(r),steps(m*b/w),N(m*b),route(N),groups(N) {
        assert(w%4==0 && m%w==0 && b%2==0 && m<=65535 && w<=65535);
        std::mt19937_64 rng(seed);std::vector<std::uint16_t> perm(m),small(w);
        for(int region=0;region<b;++region) {
            std::iota(perm.begin(),perm.end(),0);std::shuffle(perm.begin(),perm.end(),rng);
            std::copy(perm.begin(),perm.end(),route.begin()+region*m);
        }
        for(int t=0;t<steps;++t) {
            std::iota(small.begin(),small.end(),0);std::shuffle(small.begin(),small.end(),rng);
            std::copy(small.begin(),small.end(),groups.begin()+t*w);
        }
    }
};
struct Counts { i64 A=0,C=0,D=0,F=0; };
int singly_covered(const Geometry &g,const std::vector<unsigned char> &state) {
    int result=0;
    for(int t=0;t<g.steps;++t)for(int j=0;j<g.width;j+=4) {
        int offset=t*g.width,a=0;
        for(int k=0;k<4;++k)if(t)a+=state[offset-g.width+g.groups[offset+j+k]];
        for(int k=0;k<4;++k) {
            int p=offset+g.groups[offset+j+k];
            if(g.route[p]<g.inactive && a+state[p]==1)++result;
        }
    }
    return result;
}
void canonical_quartet_checks(const Geometry &g,const std::vector<unsigned char> &state,int num,int den,int carry) {
    for(int t=1;t<g.steps;++t)for(int j=0;j<g.width;j+=4) {
        int offset=t*g.width,a=0,h=0,freed=0;
        for(int k=0;k<4;++k) {
            int p=offset+g.groups[offset+j+k];
            a+=state[p-g.width];h+=state[p-g.width] && state[p];
            freed+=g.route[p]<g.inactive && !state[p];
        }
        if(a)assert(den*freed<num*a+carry*h);
    }
}
// A constructive feasible flow: serve new inactive-position arrivals first,
// then water-fill old queues within each quartet. Carry capacity is explicit.
// Lost work is an affine-expansion offset, not encoder data or code state.
i64 online_flow(const Geometry &g,int num,int den,int carry) {
    std::vector<int> previous(g.width,0),next(g.width);
    i64 served=0,lost=0;
    for(int t=0;t<g.steps;++t) {
        const int offset=t*g.width;
        for(int j=0;j<g.width;++j) {
            const bool inactive=g.route[offset+j]<g.inactive;
            next[j]=inactive?0:num;if(inactive)served+=num;
        }
        if(t)for(int j=0;j<g.width;j+=4) {
            int index[4],q[4],budget=0,total=0,maximum=0;
            for(int k=0;k<4;++k) {
                index[k]=g.groups[offset+j+k];q[k]=previous[index[k]];
                total+=q[k];maximum=std::max(maximum,q[k]);
                if(g.route[offset+index[k]]<g.inactive)budget+=den-num;
            }
            const int used=std::min(total,budget);served+=used;
            int lo=0,hi=maximum;
            while(lo<hi) {
                int middle=(lo+hi)/2,cost=0;
                for(int k=0;k<4;++k)cost+=std::max(0,q[k]-middle);
                if(cost<=used)hi=middle;else lo=middle+1;
            }
            int remaining=used;
            for(int k=0;k<4;++k){remaining-=std::max(0,q[k]-lo);q[k]=std::min(q[k],lo);}
            for(int k=0;k<4 && remaining;++k)if(q[k]==lo && q[k]){--q[k];--remaining;}
            assert(remaining==0);
            for(int k=0;k<4;++k) {
                lost+=std::max(0,q[k]-carry);next[index[k]]+=std::min(q[k],carry);
            }
        }
        previous.swap(next);
    }
    lost+=std::accumulate(previous.begin(),previous.end(),i64(0));
    assert(served+lost==i64(num)*g.N);
    return served;
}
Counts count(const Geometry &g,const std::vector<unsigned char> &state,std::vector<int> *coverage_out=nullptr) {
    Counts result;std::vector<int> row_coverage(g.rows,0);
    for(int t=0;t<g.steps;++t) {
        int offset=t*g.width;
        for(int j=0;j<g.width;++j) {
            result.A+=state[offset+j];
            if(t && state[offset-g.width+j] && !state[offset+j])++result.D;
        }
        for(int j=0;j<g.width;j+=4) {
            bool live=false;
            for(int k=0;k<4;++k)if(t && state[offset-g.width+g.groups[offset+j+k]])live=true;
            for(int k=0;k<4;++k) {
                int position=offset+g.groups[offset+j+k];
                if(live || state[position]) {
                    int row=g.route[position];++row_coverage[row];
                    if(row<g.inactive)++result.C;
                }
            }
        }
    }
    for(int h:row_coverage)result.F+=std::min(h,g.regions/2);
    if(coverage_out)*coverage_out=std::move(row_coverage);
    return result;
}

// Minimum deterministic defect charge over all full/empty temporal templates.
// Numerator units are den; no route or partition distribution is assumed.
i64 template_charge(const Geometry &g,const std::vector<unsigned char> &state,int num,int den,int carry) {
    assert(den%5==0);
    i64 previous[2]={0,std::numeric_limits<i64>::max()/4};int old_weight=0;
    i64 deaths=0;
    for(int t=0;t<g.steps;++t) {
        int weight=0;
        for(int j=0;j<g.width;++j) {
            weight+=state[t*g.width+j];
            if(t)deaths+=state[(t-1)*g.width+j] && !state[t*g.width+j];
        }
        i64 next0=i64(num)*weight+std::min(previous[0],previous[1]+i64(carry)*g.width
                                                       +i64(den/5)*(g.width-old_weight));
        i64 next1=i64(den-num)*(g.width-weight)+std::min(previous[1],previous[0]
                                                       +i64(4*den/5)*(g.width-weight));
        previous[0]=next0;previous[1]=next1;old_weight=weight;
    }
    return std::min(previous[0],previous[1])-i64(carry)*deaths;
}

// Maximal allowed state set with C_inactive=0. A state variable touches
// its own input position and the next update's quartet input positions.
std::vector<unsigned char> zero_inactive_support(const Geometry &g) {
    std::vector<unsigned char> state(g.N);
    for(int p=0;p<g.N;++p)state[p]=g.route[p]>=g.inactive;
    for(int t=1;t<g.steps;++t) {
        const int offset=t*g.width;
        for(int j=0;j<g.width;j+=4) {
            bool blocked=false;
            for(int k=0;k<4;++k)blocked|=g.route[offset+g.groups[offset+j+k]]<g.inactive;
            if(blocked)for(int k=0;k<4;++k)state[offset-g.width+g.groups[offset+j+k]]=0;
        }
    }
    return state;
}

void component_audit(const Geometry &g,const std::vector<unsigned char> &state,int num,int den,int carry,i64 expected) {
    std::vector<int> parent(g.N,-1),size(g.N,0),first(g.N,g.steps),last(g.N,-1);
    for(int p=0;p<g.N;++p)if(state[p]){parent[p]=p;size[p]=1;}
    auto root=[&](int p){while(parent[p]!=p){parent[p]=parent[parent[p]];p=parent[p];}return p;};
    auto join=[&](int a,int b){a=root(a);b=root(b);if(a==b)return;if(size[a]<size[b])std::swap(a,b);parent[b]=a;size[a]+=size[b];};
    for(int t=0;t<g.steps;++t) {
        int offset=t*g.width;
        if(t)for(int j=0;j<g.width;++j)if(state[offset+j]&&state[offset-g.width+j])join(offset+j,offset-g.width+j);
        for(int j=0;j<g.width;j+=4) {
            int incoming=-1;bool inactive=false;
            for(int k=0;k<4;++k)inactive|=g.route[offset+g.groups[offset+j+k]]<g.inactive;
            if(!inactive)continue;
            if(t)for(int k=0;k<4;++k){int p=offset-g.width+g.groups[offset+j+k];if(state[p]){if(incoming<0)incoming=p;else join(incoming,p);}}
            if(incoming>=0)for(int k=0;k<4;++k){int p=offset+g.groups[offset+j+k];if(state[p]&&g.route[p]<g.inactive)join(incoming,p);}
        }
    }
    std::vector<i64> score(g.N,0);std::vector<int> coverage(g.N,0),deaths(g.N,0);
    for(int p=0;p<g.N;++p)if(state[p]) {
        int r=root(p),t=p/g.width;score[r]+=num;first[r]=std::min(first[r],t);last[r]=std::max(last[r],t);
        if(t+1<g.steps&&!state[p+g.width]){score[r]-=carry;++deaths[r];}
    }
    for(int t=0;t<g.steps;++t)for(int j=0;j<g.width;j+=4) {
        int offset=t*g.width,incoming=-1;
        if(t)for(int k=0;k<4;++k){int p=offset-g.width+g.groups[offset+j+k];if(state[p])incoming=p;}
        for(int k=0;k<4;++k) {
            int p=offset+g.groups[offset+j+k];if(g.route[p]>=g.inactive)continue;
            int witness=state[p]?p:incoming;if(witness<0)continue;
            int r=root(witness);score[r]-=den;++coverage[r];
        }
    }
    int components=0,largest=0,largest_span=0; i64 total=0,big_score=0;int big_C=0,big_D=0;
    int paid_components=0,paid_largest=0;i64 paid_A=0,paid_C=0,paid_D=0,paid_score=0;
    std::array<i64,6> bins{},mass{},scores{};
    for(int p=0;p<g.N;++p)if(state[p]&&parent[p]==p) {
        assert(score[p]>0);++components;total+=score[p];
        if(coverage[p]){++paid_components;paid_largest=std::max(paid_largest,size[p]);paid_A+=size[p];paid_C+=coverage[p];paid_D+=deaths[p];paid_score+=score[p];}
        if(size[p]>largest){largest=size[p];largest_span=last[p]-first[p]+1;big_score=score[p];big_C=coverage[p];big_D=deaths[p];}
        int bin=size[p]<=4?0:size[p]<=16?1:size[p]<=64?2:size[p]<=256?3:size[p]<=1024?4:5;
        ++bins[bin];mass[bin]+=size[p];scores[bin]+=score[p];
    }
    assert(total==expected);
    std::cout<<"positive_components="<<components<<" largest_component_A="<<largest
             <<" largest_component_span="<<largest_span<<" largest_component_score_numerator="<<big_score
             <<" largest_component_C="<<big_C<<" largest_component_D="<<big_D<<" denominator="<<den<<std::endl;
    for(int i=0;i<6;++i)std::cout<<"component_size_bin="<<i<<" count="<<bins[i]<<" A="<<mass[i]<<" score_numerator="<<scores[i]<<std::endl;
    std::cout<<"paid_components="<<paid_components<<" paid_largest_A="<<paid_largest<<" paid_A="<<paid_A
             <<" paid_C="<<paid_C<<" paid_D="<<paid_D<<" paid_score_numerator="<<paid_score<<std::endl;
}

void hole_audit(const Geometry &g,const std::vector<unsigned char> &hole,int num,int den,int carry,i64 expected) {
    i64 birth_gate_checks=0;
    for(int t=1;t<g.steps;++t) {
        const int offset=t*g.width;
        if(t+1==g.steps){for(int j=0;j<g.width;++j)assert(!hole[offset+j]||hole[offset-g.width+j]);continue;}
        const int next=offset+g.width;
        for(int j=0;j<g.width;j+=4) {
            bool birth=false,all=true,completed=false;
            for(int k=0;k<4;++k){int p=offset+g.groups[next+j+k];birth|=hole[p]&&!hole[p-g.width];all&=hole[p];completed|=hole[p+g.width]&&g.route[p+g.width]<g.inactive;}
            if(birth){assert(all&&completed);++birth_gate_checks;}
        }
    }
    std::vector<int> parent(g.N,-1),size(g.N,0),first(g.N,g.steps),last(g.N,-1);
    for(int p=0;p<g.N;++p)if(hole[p]){parent[p]=p;size[p]=1;}
    auto root=[&](int p){while(parent[p]!=p){parent[p]=parent[parent[p]];p=parent[p];}return p;};
    auto join=[&](int a,int b){a=root(a);b=root(b);if(a==b)return;if(size[a]<size[b])std::swap(a,b);parent[b]=a;size[a]+=size[b];};
    for(int t=0;t<g.steps;++t) {
        int offset=t*g.width;
        if(t)for(int j=0;j<g.width;++j)if(hole[offset+j]&&hole[offset-g.width+j])join(offset+j,offset-g.width+j);
        if(!t)continue;
        for(int j=0;j<g.width;j+=4) {
            bool all=true;for(int k=0;k<4;++k)all&=hole[offset-g.width+g.groups[offset+j+k]];
            if(!all)continue;
            for(int k=0;k<4;++k){int p=offset+g.groups[offset+j+k];if(hole[p]&&g.route[p]<g.inactive)
                for(int l=0;l<4;++l)join(p,offset-g.width+g.groups[offset+j+l]);}
        }
    }
    std::vector<i64> score(g.N,0);std::vector<int> uncovered(g.N,0),births(g.N,0);
    for(int p=0;p<g.N;++p)if(hole[p]) {
        int r=root(p),t=p/g.width;score[r]-=num;first[r]=std::min(first[r],t);last[r]=std::max(last[r],t);
        if(t&&!hole[p-g.width]){score[r]-=carry;++births[r];}
    }
    for(int t=0;t<g.steps;++t)for(int j=0;j<g.width;j+=4) {
        int offset=t*g.width;bool all=true;
        if(t)for(int k=0;k<4;++k)all&=hole[offset-g.width+g.groups[offset+j+k]];
        if(!all)continue;
        for(int k=0;k<4;++k){int p=offset+g.groups[offset+j+k];if(hole[p]&&g.route[p]<g.inactive){int r=root(p);score[r]+=den;++uncovered[r];}}
    }
    int components=0,largest=0,span=0;i64 total=0,total_A=0,total_U=0,total_B=0,largest_score=0;
    std::array<i64,6> bins{},mass{},scores{};
    for(int p=0;p<g.N;++p)if(hole[p]&&parent[p]==p) {
        assert(score[p]>0);++components;total+=score[p];total_A+=size[p];total_U+=uncovered[p];total_B+=births[p];
        if(size[p]>largest){largest=size[p];span=last[p]-first[p]+1;largest_score=score[p];}
        int bin=size[p]<=4?0:size[p]<=16?1:size[p]<=64?2:size[p]<=256?3:size[p]<=1024?4:5;
        ++bins[bin];mass[bin]+=size[p];scores[bin]+=score[p];
    }
    assert(total==expected);
    std::cout<<"hole_components="<<components<<" hole_A="<<total_A<<" hole_U="<<total_U<<" hole_births="<<total_B
             <<" hole_gain_numerator="<<total<<" largest_hole_component_A="<<largest<<" largest_hole_span="<<span
             <<" largest_hole_gain_numerator="<<largest_score<<std::endl;
    std::cout<<"hole_birth_gate_checks="<<birth_gate_checks<<" terminal_hole_births=0"<<std::endl;
    for(int i=0;i<6;++i)std::cout<<"hole_size_bin="<<i<<" count="<<bins[i]<<" A="<<mass[i]<<" gain_numerator="<<scores[i]<<std::endl;
}

Counts run(const Geometry &g,int num,int den,bool verbose,int carry=-1,std::vector<int> *coverage_out=nullptr) {
    if(carry<0)carry=den;
    assert(carry>0);
    assert(0<num && num<den && i64(num)*g.N+1<std::numeric_limits<int>::max());
    const int source=g.N+g.N/4+g.inactive*g.regions,sink=source+1;
    const int infinity=num*g.N+1;
    const std::size_t reserve=std::size_t(8)*g.N;
    // Stored graph, original capacities, auxiliary arrays, verification
    // balances, route, groups, cuts, and component audit fit below512MiB.
    assert(reserve*sizeof(Edge)+reserve/2*sizeof(int)+std::size_t(source+2)*32
           +std::size_t(g.N)*41 < (std::size_t(512)<<20));
    Flow flow(sink+1,reserve);int next_group=g.N,next_label=g.N+g.N/4;
    for(int p=0;p<g.N;++p)flow.add(source,p,num);
    std::vector<int> label(g.width,-1);
    for(int t=0;t<g.steps;++t) {
        int offset=t*g.width;
        for(int j=0;j<g.width;++j) {
            if(t)flow.add(offset-g.width+j,offset+j,carry);
            label[j]=-1;
            if(g.route[offset+j]<g.inactive) {
                label[j]=next_label++;flow.add(offset+j,label[j],infinity);
                flow.add(label[j],sink,den);
            }
        }
        for(int j=0;j<g.width;j+=4) {
            bool needed=false;
            for(int k=0;k<4;++k)needed|=label[g.groups[offset+j+k]]>=0;
            int group=next_group++;
            if(!needed || !t)continue;
            for(int k=0;k<4;++k) {
                int coordinate=g.groups[offset+j+k];
                flow.add(offset-g.width+coordinate,group,infinity);
                if(label[coordinate]>=0)flow.add(group,label[coordinate],infinity);
            }
        }
    }
    assert(next_group==g.N+g.N/4 && next_label==source);
    if(verbose)std::cout<<"nodes="<<sink+1<<" directed_edges="<<flow.edges.size()/2<<std::endl;
    i64 value=flow.solve(source,sink,verbose);flow.verify(source,sink,value);
    const i64 constructive=online_flow(g,num,den,carry);assert(constructive<=value);
    std::vector<unsigned char> state(g.N);
    for(int p=0;p<g.N;++p)state[p]=flow.level[p]>=0;
    Counts c=count(g,state,coverage_out);
    // The residual source side is the inclusion-minimal minimum cut. If
    // an inactive input had one selected neighbor, deleting that neighbor
    // would improve the objective or preserve it with fewer state vertices.
    if(num+carry<=den)assert(singly_covered(g,state)==0);
    canonical_quartet_checks(g,state,num,den,carry);
    const i64 objective=den*c.C+carry*c.D-num*c.A;
    assert(objective==value-i64(num)*g.N);
    int internal_runs=0,terminal_runs=0,min_internal=g.steps+1;
    for(int j=0;j<g.width;++j) {
        int length=0;
        for(int t=0;t<g.steps;++t) {
            if(state[t*g.width+j])++length;
            else if(length) {
                assert(num*length>=carry);++internal_runs;
                min_internal=std::min(min_internal,length);length=0;
            }
        }
        terminal_runs+=length>0;
    }
    if(verbose) {
        component_audit(g,state,num,den,carry,-objective);
        // Reverse residual reachability selects the inclusion-minimal holes,
        // equivalently the inclusion-maximal optimal source-side support.
        std::fill(flow.level.begin(),flow.level.end(),-1);int begin=0,end=0;
        flow.queue[end++]=sink;flow.level[sink]=0;
        while(begin<end){int v=flow.queue[begin++];for(int e=flow.head[v];e!=-1;e=flow.edges[e].next)
            if(flow.edges[e^1].cap>0&&flow.level[flow.edges[e].to]<0){int w=flow.edges[e].to;flow.level[w]=0;flow.queue[end++]=w;}}
        std::vector<unsigned char> holes(g.N),maximum_state(g.N);
        for(int p=0;p<g.N;++p){holes[p]=flow.level[p]>=0;maximum_state[p]=!holes[p];}
        const auto maximal=count(g,maximum_state);
        assert(den*maximal.C+carry*maximal.D-num*maximal.A==objective);
        const i64 full_score=i64(num)*g.N-i64(den)*g.inactive*g.regions;
        hole_audit(g,holes,num,den,carry,-objective-full_score);
        const i64 W=i64(9)*g.N/20;
        std::cout<<"A="<<c.A<<" C="<<c.C<<" D="<<c.D<<" row_rank_cost="<<c.F
                 <<" objective="<<objective<<" maxflow="<<value
                 <<" affine_expansion_verified=1 slack_numerator="<<i64(num)*g.N-value
                 <<" slack_denominator="<<den
                 <<" constructive_flow_slack_numerator="<<i64(num)*g.N-constructive
                 <<" zero_slack_expansion_verified="<<(value==i64(num)*g.N)
                 <<" nonterminal_runs="<<internal_runs<<" shortest_nonterminal_run="<<min_internal
                 <<" terminal_runs="<<terminal_runs
                 <<" structural_bad_code_witness="<<(c.A>c.F && c.F+c.D<W)
                 <<" whole_ensemble_certificate=0"<<std::endl;
        if(num+carry<=den)std::cout<<"canonical_no_single_incidence_verified=1"<<std::endl;
        std::cout<<"canonical_quartet_erasure_verified=1"<<std::endl;
        if(den%5==0)std::cout<<"minimum_template_charge_numerator="<<template_charge(g,state,num,den,carry)
                            <<" denominator="<<den<<" deterministic_defect_bound=1"<<std::endl;
    }
    return c;
}

void tests() {
    Geometry g(4,8,2,2,92741);int checks=0,erasure_checks=0;
    const auto free=zero_inactive_support(g);const auto free_counts=count(g,free);
    assert(free_counts.C==0);
    for(int carry:{3,4})for(int num:{1,2,3}) {
        i64 best=0;std::vector<unsigned char> state(g.N);
        for(unsigned mask=0;mask<(1u<<g.N);++mask) {
            for(int j=0;j<g.N;++j)state[j]=(mask>>j)&1;
            auto c=count(g,state);best=std::min(best,4*c.C+carry*c.D-num*c.A);++checks;
            assert(c.A<=free_counts.A+5*c.C);
            if(c.C==0)for(int j=0;j<g.N;++j)assert(!state[j] || free[j]);
            if(mask<4096)for(int j=0;j<g.width;++j) {
                int length=0;
                for(int t=0;t<g.steps;++t) {
                    if(state[t*g.width+j])++length;
                    else if(length) {
                        if(num*length<carry) {
                            auto erased=state;
                            for(int h=t-length;h<t;++h)erased[h*g.width+j]=0;
                            auto d=count(g,erased);
                            assert(d.A==c.A-length && d.D==c.D-1 && d.C<=c.C);
                            assert(4*d.C+carry*d.D-num*d.A<4*c.C+carry*c.D-num*c.A);
                            ++erasure_checks;
                        }
                        length=0;
                    }
                }
            }
        }
        auto got=run(g,num,4,false,carry);
        assert(best==4*got.C+carry*got.D-num*got.A);
    }
    std::cout<<"exhaustive_cut_checks="<<checks<<" feasible_flow_checks=6 run_erasure_checks="<<erasure_checks
             <<" zero_neighbor_activity_checks="<<checks<<std::endl;
    // Exercise the strict and equality cases of the canonical inequality.
    for(int num:{1,2}) {
        i64 best=0;int best_size=0;
        std::vector<unsigned char> state(g.N);
        for(unsigned mask=0;mask<(1u<<g.N);++mask) {
            for(int j=0;j<g.N;++j)state[j]=(mask>>j)&1;
            auto c=count(g,state);i64 value=4*c.C+2*c.D-num*c.A;
            if(value<best || (value==best && c.A<best_size)){best=value;best_size=int(c.A);}
        }
        for(unsigned mask=0;mask<(1u<<g.N);++mask) {
            for(int j=0;j<g.N;++j)state[j]=(mask>>j)&1;
            auto c=count(g,state);
            if(4*c.C+2*c.D-num*c.A==best && c.A==best_size) {
                assert(singly_covered(g,state)==0);
                canonical_quartet_checks(g,state,num,4,2);
            }
        }
        run(g,num,4,false,2);
    }
    std::cout<<"canonical_toy_cases=2 strict_and_tie_cases=1"<<std::endl;
    // All supports and all temporal templates on the small geometry.
    std::vector<Counts> baseline(1u<<g.steps);
    for(unsigned bits=0;bits<baseline.size();++bits) {
        std::vector<unsigned char> state(g.N);
        for(int t=0;t<g.steps;++t)for(int j=0;j<g.width;++j)state[t*g.width+j]=(bits>>t)&1;
        baseline[bits]=count(g,state);
    }
    i64 defect_checks=0;
    for(unsigned mask=0;mask<(1u<<g.N);++mask) {
        std::vector<unsigned char> state(g.N);
        for(int j=0;j<g.N;++j)state[j]=(mask>>j)&1;
        auto c=count(g,state);i64 best=std::numeric_limits<i64>::max();
        for(unsigned bits=0;bits<baseline.size();++bits) {
            int plus=0,minus=0,starts=0,ends=0;
            for(int t=0;t<g.steps;++t) {
                int weight=std::popcount((mask>>(t*g.width))&((1u<<g.width)-1));
                if((bits>>t)&1) {
                    minus+=g.width-weight;
                    if(!t || !((bits>>(t-1))&1))starts+=g.width-weight;
                    if(t+1<g.steps && !((bits>>(t+1))&1))ends+=g.width-weight;
                } else plus+=weight;
            }
            auto b=baseline[bits];
            assert(5*(b.C-c.C)<=5*minus+4*starts+ends);
            i64 charge=25*plus+75*minus+80*starts+20*ends-74*(c.D-b.D);
            assert(25*c.A-100*c.C-74*c.D<=25*b.A-100*b.C-74*b.D+charge);
            best=std::min(best,charge);++defect_checks;
        }
        assert(best==template_charge(g,state,25,100,74));
    }
    std::cout<<"exact_template_defect_checks="<<defect_checks<<" exact_template_dp_checks="<<(1u<<g.N)<<std::endl;
}

int main(int argc,char **argv) {
    const auto start=std::chrono::steady_clock::now();tests();
    if(argc>1 && std::string(argv[1])=="--tests-only")return 0;
    int num=argc>1?std::stoi(argv[1]):22,den=100;
    std::uint64_t seed=argc>2?std::stoull(argv[2]):20261027;
    int carry=argc>3?std::stoi(argv[3]):den;
    Geometry g(512,8192,256,2100,seed);
    const auto free_counts=count(g,zero_inactive_support(g));assert(free_counts.C==0);
    std::cout<<"maximal_zero_inactive_activity="<<free_counts.A<<" deaths="<<free_counts.D
             <<" row_rank_cost="<<free_counts.F<<" all_support_bound=A_le_Z_plus_5C"
             <<" ensemble_tail_proved=0"<<std::endl;
    std::cout<<"K=1048576 rate=0.5 target_distance=0.45 inactive_rows=2100 seed="<<seed
             <<" slope="<<num<<"/"<<den<<" death_weight="<<carry<<"/"<<den
             <<" pseudorandom_geometry_diagnostic=1"<<std::endl;
    const int rounds=argc>4?std::stoi(argv[4]):1;assert(rounds>=1&&rounds<=3);
    for(int round=0;round<rounds;++round) {
        std::vector<int> coverage;const auto c=run(g,num,den,true,carry,&coverage);
        std::vector<int> order(g.rows),relabel(g.rows);std::iota(order.begin(),order.end(),0);
        std::sort(order.begin(),order.end(),[&](int a,int b){return coverage[a]!=coverage[b]?coverage[a]<coverage[b]:a<b;});
        i64 optimized_C=0;int new_rows=0;
        for(int j=0;j<g.inactive;++j){optimized_C+=coverage[order[j]];new_rows+=order[j]>=g.inactive;}
        assert(optimized_C<=c.C);
        std::cout<<"row_response_round="<<round<<" A="<<c.A<<" D="<<c.D<<" old_C="<<c.C<<" chosen_C="<<optimized_C
                 <<" new_rows="<<new_rows<<" score_lower_numerator="<<i64(num)*c.A-i64(carry)*c.D-i64(den)*optimized_C
                 <<" denominator="<<den<<" meets_residual_activity="<<(c.A>=743225)
                 <<" meets_residual_turnover="<<(optimized_C+c.D<=255463)
                 <<" meets_residual_rank_score="<<(c.A-optimized_C>=713983)
                 <<" geometry_only=1"<<std::endl;
        if(optimized_C==c.C||round+1==rounds)break;
        for(int j=0;j<g.rows;++j)relabel[order[j]]=j;
        for(auto &row:g.route)row=static_cast<std::uint16_t>(relabel[row]);
    }
    std::cout<<"seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
             <<" files_written=0"<<std::endl;
    return 0;
}
