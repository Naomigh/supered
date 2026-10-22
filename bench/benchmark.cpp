#include "../src/references.hpp"
#include "strong/myers_v1.hpp"
#include "strong/myers_v2.hpp"
#include "strong/banded.hpp"
#include "supered/supered.hpp"
#include "exploration.hpp"
#include "perf_gate.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <sched.h>
#include <sys/resource.h>
#include <time.h>

using namespace strong;
using namespace test_support;
namespace fs=std::filesystem;
constexpr int base_refs=4, max_refs=20;
const char* algorithm_names[]={"edlib_original_hot","edlib_original_stack",
    "myers_banded64x4","myers_128x4","leap_corrected_avx2","myers_2x64x4_v1",
    "myers_2x64x4_v2","banded_u1","banded_u2","banded_u4","banded_u8","supered"};
struct alignas(32) RawTile {
    PackedDNA128 fixed[4],stream[max_refs];
    int16_t exact[4*max_refs];
};
struct Scratch {
    Sequence fixed[4],stream[max_refs];
    void load(const RawTile& raw,int nq,int nt) {
        for(int i=0;i<nq;++i) fixed[i]=Sequence(decode(raw.fixed[i]));
        for(int i=0;i<nt;++i) stream[i]=Sequence(decode(raw.stream[i]));
    }
};
inline uint64_t cpu_ns() {
    timespec t;clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t);
    return uint64_t(t.tv_sec)*1000000000+t.tv_nsec;
}
inline uint64_t wall_ns() {
    timespec t;clock_gettime(CLOCK_MONOTONIC_RAW,&t);
    return uint64_t(t.tv_sec)*1000000000+t.tv_nsec;
}
inline void touch(const void* ptr,size_t bytes) {
    auto* p=(const volatile unsigned char*)ptr;
    unsigned sum=0;
    for(size_t i=0;i<bytes;i+=64) sum+=p[i];
    asm volatile(""::"r"(sum):"memory");_mm_mfence();
}
template<int K,int Algorithm>
__attribute__((noinline)) uint64_t compute(const Scratch& s,int nq,int nt,SIMD_ED& leap,
                                         supered::Context<K>& combined) {
    // One call per timing window. Pair kernels remain inline/optimizable.
    asm volatile("":::"memory");
    uint64_t sum=0;
    for(int q=0;q<nq;++q) for(int t=0;t<nt;t+=4) {
        if constexpr(Algorithm==2 || Algorithm==3 || Algorithm>=5) {
            Quad query=broadcast(s.fixed[q].packed);
            Quad target{&s.stream[t].packed,&s.stream[t+1].packed,
                        &s.stream[t+2].packed,&s.stream[t+3].packed};
            Distances d;
            if constexpr(Algorithm==2) d=myers_banded64x4<K>(query,target);
            else if constexpr(Algorithm==3) d=myers_128x4<K>(query,target);
            else if constexpr(Algorithm==5) d=v1::myers_2x64x4<K>(query,target);
            else if constexpr(Algorithm==6) d=v2::myers_2x64x4<K>(query,target);
            else if constexpr(Algorithm==7) d=myers_banded64x4<K,1>(query,target);
            else if constexpr(Algorithm==8) d=myers_banded64x4<K,2>(query,target);
            else if constexpr(Algorithm==9) d=myers_banded64x4<K,4>(query,target);
            else if constexpr(Algorithm==11) d=combined.distance4(query,target);
            else d=myers_banded64x4<K,8>(query,target);
            for(int j=0;j<4;++j) sum+=uint64_t(d[j]+1);
        } else for(int j=0;j<4;++j) {
            int d;
            if constexpr(Algorithm==0) d=reference::hot<K,false>(s.fixed[q],s.stream[t+j]);
            else if constexpr(Algorithm==1) d=reference::hot<K,true>(s.fixed[q],s.stream[t+j]);
            else d=reference::leap(leap,s.fixed[q].packed,s.stream[t+j].packed);
            sum+=uint64_t(d+1);
        }
    }
    asm volatile(""::"r"(sum):"memory");return sum;
}
template<int K>
void verify_tile(const Scratch& s,int nq,RawTile& raw,SIMD_ED& leap,
                 std::map<int,uint64_t>& histogram,supered::Context<K>& combined) {
    for(int q=0;q<nq;++q) for(int t=0;t<max_refs;t+=4) {
        Quad query=broadcast(s.fixed[q].packed);
        Quad target{&s.stream[t].packed,&s.stream[t+1].packed,
                    &s.stream[t+2].packed,&s.stream[t+3].packed};
        auto a=myers_banded64x4<K>(query,target),b=myers_128x4<K>(query,target);
        auto routed=combined.distance4(query,target);
        auto fa=myers_2x64x4<128>(query,target),fb=myers_128x4<128>(query,target);
        for(int j=0;j<4;++j) {
            int d=reference::api(s.fixed[q],s.stream[t+j]);
            raw.exact[q*max_refs+t+j]=d;++histogram[d];
            int expected=d<=K?d:-1;
            int old=reference::hot<K,false>(s.fixed[q],s.stream[t+j]);
            int stack=reference::hot<K,true>(s.fixed[q],s.stream[t+j]);
            int lp=K<=5?reference::leap(leap,s.fixed[q].packed,s.stream[t+j].packed):expected;
            if(routed[j]!=expected||a[j]!=expected||b[j]!=expected||old!=expected||stack!=expected||lp!=expected
                ||fa[j]!=d||fb[j]!=d) throw std::runtime_error("benchmark pair oracle mismatch");
        }
    }
}
struct Totals {
    uint64_t cpu=0,wall=0,cycles=0,l2=0,sum=0,rejected=0,reject_ns=0,reject_dram=0,max_cpu=0;
};
struct WindowRecord {
    uint64_t cpu_ns,wall_ns,user_cycles;
    uint32_t l2_misses,attempts;
};
static_assert(sizeof(WindowRecord)==32);
template<int K,int Algorithm>
void trial(const std::vector<RawTile>& data,int nq,int multiplier,int repeat,
           const fs::path& out,std::ofstream& trials,std::ofstream& rejected) {
    int nt=base_refs*multiplier;
    Scratch scratch;SIMD_ED leap;supered::Context<K> combined;
    if constexpr(Algorithm==4) leap.init_levenshtein(K,ED_GLOBAL,false);
    scratch.load(data[0],nq,nt);
    compute<K,Algorithm>(scratch,nq,nt,leap,combined); // code/state warmup, outside timers
    std::vector<WindowRecord> windows;windows.reserve(data.size());
    PerfGate pmu;
    rusage before{},after{};getrusage(RUSAGE_THREAD,&before);
    uint64_t expected_sum=0,passed=0;Totals totals;
    pmu.start();
    for(size_t w=0;w<data.size();++w) {
        const auto& raw=data[w];
        scratch.load(raw,nq,nt); // stream loading, decoding and profiles EXCLUDED
        uint64_t expected=0;
        for(int q=0;q<nq;++q) for(int t=0;t<nt;++t) {
            int d=raw.exact[q*max_refs+t];
            if(d<=K) {expected+=d+1;++passed;}
        }
        expected_sum+=expected;bool accepted=false;
        for(int attempt=0;attempt<1024;++attempt) {
            touch(scratch.fixed,size_t(nq)*sizeof(Sequence));
            touch(scratch.stream,size_t(nt)*sizeof(Sequence));
            auto w0=wall_ns();auto t0=cpu_ns();auto p0=pmu.snapshot();
            uint64_t sum=compute<K,Algorithm>(scratch,nq,nt,leap,combined);
            auto p1=pmu.snapshot();auto t1=cpu_ns();auto w1=wall_ns();
            uint64_t cpu=t1-t0,wall=w1-w0,cycles=p1[0]-p0[0],l2=p1[1]-p0[1],dram=p1[2]-p0[2];
            if(sum!=expected) throw std::runtime_error("timed checksum mismatch");
            if(dram) {
                ++totals.rejected;totals.reject_ns+=cpu;totals.reject_dram+=dram;
                rejected<<K<<','<<algorithm_names[Algorithm]<<','<<repeat<<','<<multiplier<<','
                    <<w<<','<<attempt<<','<<cpu<<','<<cycles<<','<<l2<<','<<dram<<'\n';
                continue;
            }
            totals.cpu+=cpu;totals.wall+=wall;totals.cycles+=cycles;totals.l2+=l2;totals.sum+=sum;
            totals.max_cpu=std::max(totals.max_cpu,cpu);
            windows.push_back({cpu,wall,cycles,uint32_t(l2),uint32_t(attempt+1)});
            accepted=true;break;
        }
        if(!accepted) throw std::runtime_error("cache gate exhausted: input was not skipped");
    }
    pmu.stop();getrusage(RUSAGE_THREAD,&after);
    if(totals.sum!=expected_sum) throw std::runtime_error("aggregate checksum mismatch");
    std::string stem=std::string(algorithm_names[Algorithm])+"_r"+std::to_string(repeat)+"_m"+std::to_string(multiplier);
    std::ofstream binary(out/(stem+".windows.bin"),std::ios::binary);
    binary.write((const char*)windows.data(),windows.size()*sizeof(WindowRecord));
    uint64_t n=uint64_t(data.size())*nq*nt;
    trials<<K<<','<<algorithm_names[Algorithm]<<','<<repeat<<','<<multiplier<<','<<n<<','
        <<data.size()<<','<<nq<<','<<nt<<','<<sizeof(Sequence)*(nq+nt)<<','
        <<totals.cpu<<','<<totals.wall<<','<<totals.cycles<<','<<totals.l2<<",0,"
        <<totals.rejected<<','<<totals.reject_ns<<','<<totals.reject_dram<<','
        <<totals.sum<<','<<expected_sum<<','<<passed<<','<<pmu.enabled<<','<<pmu.running<<','
        <<after.ru_minflt-before.ru_minflt<<','<<after.ru_majflt-before.ru_majflt<<','
        <<after.ru_nivcsw-before.ru_nivcsw<<','<<totals.max_cpu<<'\n';
    trials.flush();rejected.flush();
    std::cout<<"TRIAL K="<<K<<" algorithm="<<algorithm_names[Algorithm]<<" repeat="<<repeat
        <<" N="<<n<<" cpu_ns="<<totals.cpu<<" rejected="<<totals.rejected<<std::endl;
}
template<int K>
void run(const fs::path& out,const std::string& workload,size_t ntiles,int repeats,const fs::path& existing) {
    if(fs::exists(out/"trials.csv")) throw std::runtime_error("refusing to overwrite trials");
    fs::create_directories(out);
    int nq=workload=="near"?1:4;
    uint64_t seed=0x2026102200ULL+K*100+nq;
    std::mt19937_64 rng(seed);
    std::vector<RawTile> data(ntiles);Scratch scratch;SIMD_ED leap;supered::Context<K> combined;
    if(!existing.empty()) {
        std::ifstream input(existing,std::ios::binary);
        input.read((char*)data.data(),data.size()*sizeof(RawTile));
        if(!input) throw std::runtime_error("input corpus read");
    }
    if constexpr(K<=5) leap.init_levenshtein(K,ED_GLOBAL,false);
    std::map<int,uint64_t> histogram;
    for(size_t w=0;w<ntiles;++w) {
      if(!existing.empty()) scratch.load(data[w],nq,max_refs);
      else {
        for(int q=0;q<nq;++q) {
            scratch.fixed[q]=Sequence(random_dna(rng));data[w].fixed[q]=scratch.fixed[q].packed;
        }
        for(int t=0;t<max_refs;++t) {
            DNA dna;
            if(workload=="cross" && t%2) dna=random_dna(rng);
            else {
                int source=rng()%nq,edits=rng()%(K+1);
                dna=mutate(scratch.fixed[source].dna,edits,rng,(w+t)%2);
            }
            scratch.stream[t]=Sequence(dna);data[w].stream[t]=scratch.stream[t].packed;
        }
      }
        verify_tile<K>(scratch,nq,data[w],leap,histogram,combined);
        if((w+1)%8192==0) std::cout<<"ORACLE K="<<K<<" workload="<<workload<<" tiles="<<w+1<<'/'<<ntiles<<std::endl;
    }
    std::ofstream inputs(out/"inputs.bin",std::ios::binary);
    inputs.write((const char*)data.data(),data.size()*sizeof(RawTile));inputs.close();
    std::ofstream meta(out/"dataset.json");
    meta<<"{\"K\":"<<K<<",\"workload\":\""<<workload<<"\",\"seed\":"<<seed
        <<",\"tiles\":"<<ntiles<<",\"fixed_per_tile\":"<<nq<<",\"max_stream_per_tile\":20,"
        <<"\"raw_tile_bytes\":"<<sizeof(RawTile)<<",\"sequence_bytes\":"<<sizeof(Sequence)
        <<",\"base_refs\":4,\"repeats\":"<<repeats<<",\"all_max_point_pairs_validated\":"
        <<uint64_t(ntiles)*nq*max_refs<<"}\n";
    std::ofstream hist(out/"distance_histogram.csv");hist<<"distance,pairs\n";
    for(auto [d,n]:histogram) hist<<d<<','<<n<<'\n';
    std::ofstream trials(out/"trials.csv"),rejected(out/"rejected_windows.csv");
    trials<<"K,algorithm,repeat,multiplier,N,windows,fixed_sequences,stream_sequences,prepared_bytes,"
        "cpu_ns,wall_ns,user_cycles,accepted_l2_misses,accepted_dram_fills,rejected_windows,"
        "rejected_cpu_ns,rejected_dram_fills,checksum,expected_checksum,passing_pairs,enabled_ns,"
        "running_ns,minor_faults,major_faults,involuntary_switches,max_window_cpu_ns\n";
    rejected<<"K,algorithm,repeat,multiplier,window,attempt,cpu_ns,user_cycles,l2_misses,dram_fills\n";
    int algorithms=7;
    std::vector<int> selected;
    if(const char* choices=std::getenv("STRONG_ALGORITHMS")) {
        std::stringstream ss(choices);std::string s;
        while(std::getline(ss,s,',')) selected.push_back(std::stoi(s));
    } else for(int a=0;a<algorithms;++a) selected.push_back(a);
    for(int r=0;r<repeats;++r) {
        std::vector<std::pair<int,int>> jobs;
        for(int a:selected) {
            if(K>5 && a==4) continue;
            for(int m=1;m<=5;++m) jobs.emplace_back(a,m);
        }
        std::shuffle(jobs.begin(),jobs.end(),rng);
        for(auto [a,m]:jobs) switch(a) {
            case 0:trial<K,0>(data,nq,m,r,out,trials,rejected);break;
            case 1:trial<K,1>(data,nq,m,r,out,trials,rejected);break;
            case 2:trial<K,2>(data,nq,m,r,out,trials,rejected);break;
            case 3:trial<K,3>(data,nq,m,r,out,trials,rejected);break;
            case 4:if constexpr(K<=5) trial<K,4>(data,nq,m,r,out,trials,rejected);break;
            case 5:trial<K,5>(data,nq,m,r,out,trials,rejected);break;
            case 6:trial<K,6>(data,nq,m,r,out,trials,rejected);break;
            case 7:trial<K,7>(data,nq,m,r,out,trials,rejected);break;
            case 8:trial<K,8>(data,nq,m,r,out,trials,rejected);break;
            case 9:trial<K,9>(data,nq,m,r,out,trials,rejected);break;
            case 10:trial<K,10>(data,nq,m,r,out,trials,rejected);break;
            case 11:trial<K,11>(data,nq,m,r,out,trials,rejected);break;
        }
    }
    std::ofstream(out/"COMPLETE")<<"All pairs validated; all accepted windows DRAM/IO demand fills = 0.\n";
}
int main(int argc,char**argv) try {
    if(argc!=6&&argc!=7) throw std::runtime_error("OUTPUT K near|cross TILES REPEATS [EXISTING_INPUTS]");
    cpu_set_t set;CPU_ZERO(&set);CPU_SET(24,&set);
    if(sched_setaffinity(0,sizeof(set),&set)) throw std::runtime_error("affinity");
    int k=std::stoi(argv[2]),repeats=std::stoi(argv[5]);
    size_t tiles=std::stoull(argv[4]);std::string workload=argv[3];
    fs::path existing=argc==7?argv[6]:"";
    if(workload!="near"&&workload!="cross") throw std::runtime_error("workload");
    switch(k) {
        case 0:run<0>(argv[1],workload,tiles,repeats,existing);break;
        case 2:run<2>(argv[1],workload,tiles,repeats,existing);break;
        case 5:run<5>(argv[1],workload,tiles,repeats,existing);break;
        case 50:run<50>(argv[1],workload,tiles,repeats,existing);break;
        default:throw std::runtime_error("K must be 0,2,5,50");
    }
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
