// Reuse the frozen four-limb kernels without changing their hashed source.
#define main frozen_pair_search_main
#include "pdual_low_weight.cpp"
#undef main

int main(int argc,char**argv) {
    if(argc!=3)return 2;
    const double seconds=std::stod(argv[1]);const uint64_t seed=std::stoull(argv[2]);
    std::array<Word,125> original,rows;std::array<unsigned,125> tags,original_tags,projection;
    for(unsigned i=0;i<125;++i){std::string s;if(!(std::cin>>s>>original_tags[i]))return 3;original[i]=parse(s);}
    std::array<unsigned,256> columns;for(unsigned i=0;i<256;++i)columns[i]=i;
    std::array<bool,256> pivot_column;
    std::array<int,1024> head;
    std::array<int,1891> next;
    std::array<Word,1891> left;
    std::array<unsigned,1891> left_tag;
    std::array<unsigned,10> parity_columns;
    std::mt19937_64 random(seed);
    Word best_word=original[0],best_nonzero_word=original[0];
    unsigned best=256,best_nonzero=256;uint64_t trials=0,candidates=0,four_candidates=0;
    const auto start=std::chrono::steady_clock::now();
    auto elapsed=[&](){return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
    auto consider=[&](Word w,unsigned tag){
        ++candidates;unsigned count=weight(w);
        if(count>128){w=complement(w);count=256-count;}
        if(count==0)return;
        if(count<best){best=count;best_word=w;}
        if(tag && count<best_nonzero){best_nonzero=count;best_nonzero_word=w;}
    };
    do {
        rows=original;tags=original_tags;pivot_column.fill(false);
        std::shuffle(columns.begin(),columns.end(),random);
        unsigned rank=0;
        for(unsigned col:columns){
            unsigned j=rank;while(j<125 && !bit(rows[j],col))++j;
            if(j==125)continue;
            std::swap(rows[rank],rows[j]);std::swap(tags[rank],tags[j]);
            const Word pivot=rows[rank];const unsigned tag=tags[rank];pivot_column[col]=true;
            for(unsigned i=0;i<125;++i)if(i!=rank && bit(rows[i],col)){
                rows[i]=add(rows[i],pivot);tags[i]^=tag;
            }
            if(++rank==125)break;
        }
        if(rank!=125)return 4;
        unsigned selected=0;
        for(unsigned col:columns)if(!pivot_column[col]){
            parity_columns[selected++]=col;if(selected==10)break;
        }
        for(unsigned i=0;i<125;++i){
            consider(rows[i],tags[i]);
            for(unsigned j=i+1;j<125;++j)consider(add(rows[i],rows[j]),tags[i]^tags[j]);
            unsigned key=0;
            for(unsigned j=0;j<10;++j)key|=unsigned(bit(rows[i],parity_columns[j]))<<j;
            projection[i]=key;
        }
        head.fill(-1);unsigned used=0;
        for(unsigned i=0;i<62;++i)for(unsigned j=i+1;j<62;++j){
            const unsigned key=projection[i]^projection[j];
            left[used]=add(rows[i],rows[j]);left_tag[used]=tags[i]^tags[j];
            next[used]=head[key];head[key]=used++;
        }
        if(used!=1891)return 5;
        for(unsigned i=62;i<125;++i)for(unsigned j=i+1;j<125;++j){
            const unsigned key=projection[i]^projection[j];
            const Word right=add(rows[i],rows[j]);const unsigned tag=tags[i]^tags[j];
            for(int k=head[key];k!=-1;k=next[k]){
                consider(add(left[k],right),left_tag[k]^tag);++four_candidates;
            }
        }
        ++trials;
    }while(best>30 && elapsed()<seconds);
    std::cout<<"{\"trials\":"<<trials<<",\"candidates\":"<<candidates
        <<",\"four_row_candidates\":"<<four_candidates<<",\"elapsed_seconds\":"<<elapsed()
        <<",\"best_weight\":"<<best<<",\"best_word\":\""<<hex(best_word)
        <<"\",\"best_nonzero_F7_weight\":"<<best_nonzero
        <<",\"best_nonzero_F7_word\":\""<<hex(best_nonzero_word)<<"\"}\n";
}
