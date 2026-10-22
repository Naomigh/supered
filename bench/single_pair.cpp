#include "../src/references.hpp"
#include "supered/supered.hpp"
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
constexpr int base_pairs=32,max_pairs=5*base_pairs;
constexpr const char* names[]={"supered_single","supered_batch4","edlib_original_hot","edlib_stack_reference"};
struct RawPair { PackedDNA128 q,t; };
static_assert(sizeof(RawPair)==64);
struct PreparedPair { Sequence q,t; };
struct Scratch {
    Sequence q[max_pairs],t[max_pairs];
    void load(const PreparedPair* pairs,int count) {
        for(int i=0;i<count;++i) {
            q[i]=pairs[i].q;
            t[i]=pairs[i].t;
        }
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
inline void touch(const void* pointer,size_t bytes) {
    const auto* p=(const volatile unsigned char*)pointer;
    unsigned sum=0;
    for(size_t i=0;i<bytes;i+=64)sum+=p[i];
    asm volatile(""::"r"(sum):"memory");_mm_mfence();
}
template<int K,int Method>
__attribute__((noinline)) uint64_t compute(const Scratch& s,int count,supered::Context<K>& context) {
    // One boundary per timing window; every public pair kernel remains inline.
    asm volatile("":::"memory");
    uint64_t sum=0;
    if constexpr(Method==1) {
        for(int i=0;i<count;i+=4) {
            auto d=context.distance4(Quad{&s.q[i].packed,&s.q[i+1].packed,&s.q[i+2].packed,&s.q[i+3].packed},
                                     Quad{&s.t[i].packed,&s.t[i+1].packed,&s.t[i+2].packed,&s.t[i+3].packed});
            for(int lane=0;lane<4;++lane)sum+=uint64_t(d[lane]+1);
        }
    } else {
        for(int i=0;i<count;++i) {
            int d;
            if constexpr(Method==0)d=context.distance(s.q[i].packed,s.t[i].packed);
            else if constexpr(Method==2)d=reference::hot<K,false>(s.q[i],s.t[i]);
            else d=reference::hot<K,true>(s.q[i],s.t[i]);
            sum+=uint64_t(d+1);
        }
    }
    asm volatile(""::"r"(sum):"memory");return sum;
}
struct WindowRecord {uint64_t cpu,wall,cycles;uint32_t l2,attempts;};
static_assert(sizeof(WindowRecord)==32);

template<int K,int Method>
void trial(const std::vector<PreparedPair>& data,const std::vector<int16_t>& exact,size_t windows,
           int multiplier,int repeat,const fs::path& out,std::ofstream& trials,std::ofstream& rejected) {
    int count=base_pairs*multiplier;
    Scratch scratch;supered::Context<K> context;
    scratch.load(data.data(),count);
    compute<K,Method>(scratch,count,context);
    std::vector<WindowRecord> records;records.reserve(windows);
    uint64_t total_cpu=0,total_wall=0,total_cycles=0,total_l2=0,total_sum=0,expected_sum=0;
    uint64_t rejects=0,reject_ns=0,reject_dram=0,passing=0;
    PerfGate pmu;rusage before{},after{};getrusage(RUSAGE_THREAD,&before);pmu.start();
    for(size_t w=0;w<windows;++w) {
        size_t offset=w*max_pairs;
        scratch.load(data.data()+offset,count); // all preparation and cache filling excluded
        uint64_t expected=0;
        for(int i=0;i<count;++i)if(exact[offset+i]<=K){expected+=exact[offset+i]+1;++passing;}
        expected_sum+=expected;
        bool accepted=false;
        for(int attempt=0;attempt<1024;++attempt) {
            touch(scratch.q,size_t(count)*sizeof(Sequence));
            touch(scratch.t,size_t(count)*sizeof(Sequence));
            auto w0=wall_ns();auto t0=cpu_ns();auto p0=pmu.snapshot();
            uint64_t sum=compute<K,Method>(scratch,count,context);
            auto p1=pmu.snapshot();auto t1=cpu_ns();auto w1=wall_ns();
            auto cpu=t1-t0,wall=w1-w0,cycles=p1[0]-p0[0],l2=p1[1]-p0[1],dram=p1[2]-p0[2];
            if(sum!=expected)throw std::runtime_error("timed checksum mismatch");
            if(dram) {
                ++rejects;reject_ns+=cpu;reject_dram+=dram;
                rejected<<K<<','<<names[Method]<<','<<repeat<<','<<multiplier<<','<<w<<','<<attempt
                    <<','<<cpu<<','<<cycles<<','<<l2<<','<<dram<<'\n';
                continue;
            }
            total_cpu+=cpu;total_wall+=wall;total_cycles+=cycles;total_l2+=l2;total_sum+=sum;
            records.push_back({cpu,wall,cycles,uint32_t(l2),uint32_t(attempt+1)});
            accepted=true;break;
        }
        if(!accepted)throw std::runtime_error("cache gate exhausted without skipping input");
    }
    pmu.stop();getrusage(RUSAGE_THREAD,&after);
    if(total_sum!=expected_sum)throw std::runtime_error("aggregate checksum mismatch");
    std::string stem=std::string(names[Method])+"_r"+std::to_string(repeat)+"_m"+std::to_string(multiplier);
    std::ofstream binary(out/(stem+".windows.bin"),std::ios::binary);
    binary.write((const char*)records.data(),records.size()*sizeof(WindowRecord));
    trials<<K<<','<<names[Method]<<','<<repeat<<','<<multiplier<<','<<windows*count<<','<<windows
        <<','<<count<<','<<2*count*sizeof(Sequence)<<','<<total_cpu<<','<<total_wall<<','<<total_cycles
        <<','<<total_l2<<",0,"<<rejects<<','<<reject_ns<<','<<reject_dram<<','<<total_sum<<','<<expected_sum
        <<','<<passing<<','<<pmu.enabled<<','<<pmu.running<<','<<after.ru_minflt-before.ru_minflt
        <<','<<after.ru_majflt-before.ru_majflt<<','<<after.ru_nivcsw-before.ru_nivcsw<<'\n';
    trials.flush();rejected.flush();
    std::cout<<"TRIAL K="<<K<<" method="<<names[Method]<<" repeat="<<repeat<<" N="<<windows*count
        <<" cpu_ns="<<total_cpu<<" rejected="<<rejects<<std::endl;
}

template<int K>
void run(const fs::path& out,const std::string& workload,size_t windows,int repeats) {
    if(fs::exists(out/"trials.csv"))throw std::runtime_error("refusing to overwrite results");
    fs::create_directories(out);
    uint64_t seed=0x202610230000ULL+K*100+(workload=="mixed");
    std::mt19937_64 rng(seed);
    std::vector<RawPair> data(windows*max_pairs);
    std::vector<PreparedPair> prepared(data.size());
    std::vector<int16_t> exact(data.size());
    supered::Context<K> context;Scratch scratch;std::map<int,uint64_t> histogram;
    uint64_t dp_count=0,checksum=0;
    for(size_t w=0;w<windows;++w) {
        for(int i=0;i<max_pairs;++i) {
            auto a=random_dna(rng);
            auto b=workload=="mixed"&&i%2?random_dna(rng):mutate(a,rng()%(K+1),rng,(w+i)%2);
            scratch.q[i]=Sequence(a);scratch.t[i]=Sequence(b);
            data[w*max_pairs+i]={scratch.q[i].packed,scratch.t[i].packed};
            prepared[w*max_pairs+i]={scratch.q[i],scratch.t[i]};
            int d=reference::api(scratch.q[i],scratch.t[i]);
            exact[w*max_pairs+i]=d;++histogram[d];checksum+=d;
            if(i==int(w%max_pairs)) {
                if(scalar_dp(a,b)!=d)throw std::runtime_error("independent DP mismatch");
                ++dp_count;
            }
            int expected=d<=K?d:-1;
            if(context.distance(scratch.q[i].packed,scratch.t[i].packed)!=expected
                ||reference::hot<K,false>(scratch.q[i],scratch.t[i])!=expected
                ||reference::hot<K,true>(scratch.q[i],scratch.t[i])!=expected)
                throw std::runtime_error("single/original Edlib pair mismatch");
            if(workload=="near"&&d>K)throw std::runtime_error("near generator exceeded budget");
        }
        for(int i=0;i<max_pairs;i+=4) {
            auto d=context.distance4(Quad{&scratch.q[i].packed,&scratch.q[i+1].packed,&scratch.q[i+2].packed,&scratch.q[i+3].packed},
                                     Quad{&scratch.t[i].packed,&scratch.t[i+1].packed,&scratch.t[i+2].packed,&scratch.t[i+3].packed});
            for(int lane=0;lane<4;++lane) {
                int e=exact[w*max_pairs+i+lane];
                if(d[lane]!=(e<=K?e:-1))throw std::runtime_error("four-pair/original Edlib mismatch");
            }
        }
        if((w+1)%256==0)std::cout<<"ORACLE K="<<K<<" workload="<<workload<<" windows="<<w+1<<'/'<<windows<<std::endl;
    }
    std::ofstream input(out/"inputs.bin",std::ios::binary);
    input.write((const char*)data.data(),data.size()*sizeof(RawPair));input.close();
    std::ofstream distances(out/"exact_distances.i16",std::ios::binary);
    distances.write((const char*)exact.data(),exact.size()*sizeof(int16_t));distances.close();
    std::ofstream histogram_file(out/"distance_histogram.csv");histogram_file<<"distance,pairs\n";
    for(auto [d,n]:histogram)histogram_file<<d<<','<<n<<'\n';
    std::ofstream meta(out/"dataset.json");
    meta<<"{\"K\":"<<K<<",\"workload\":\""<<workload<<"\",\"seed\":"<<seed<<",\"windows\":"<<windows
        <<",\"base_pairs\":"<<base_pairs<<",\"max_pairs\":"<<max_pairs<<",\"pair_bytes\":"<<sizeof(RawPair)
        <<",\"Sequence_bytes\":"<<sizeof(Sequence)<<",\"repeats\":"<<repeats<<",\"oracle_pairs\":"<<data.size()
        <<",\"independent_DP_pairs\":"<<dp_count<<",\"distance_checksum\":"<<checksum<<"}\n";
    std::ofstream trials(out/"trials.csv"),rejected(out/"rejected_windows.csv");
    trials<<"K,algorithm,repeat,multiplier,N,windows,pairs_per_window,prepared_bytes,cpu_ns,wall_ns,user_cycles,"
        "accepted_l2_misses,accepted_dram_fills,rejected_windows,rejected_cpu_ns,rejected_dram_fills,checksum,"
        "expected_checksum,passing_pairs,enabled_ns,running_ns,minor_faults,major_faults,involuntary_switches\n";
    rejected<<"K,algorithm,repeat,multiplier,window,attempt,cpu_ns,user_cycles,l2_misses,dram_fills\n";
    for(int repeat=0;repeat<repeats;++repeat) {
        std::vector<std::pair<int,int>> jobs;
        for(int method=0;method<4;++method)for(int m=1;m<=5;++m)jobs.emplace_back(method,m);
        std::shuffle(jobs.begin(),jobs.end(),rng);
        for(auto [method,m]:jobs)switch(method) {
            case 0:trial<K,0>(prepared,exact,windows,m,repeat,out,trials,rejected);break;
            case 1:trial<K,1>(prepared,exact,windows,m,repeat,out,trials,rejected);break;
            case 2:trial<K,2>(prepared,exact,windows,m,repeat,out,trials,rejected);break;
            case 3:trial<K,3>(prepared,exact,windows,m,repeat,out,trials,rejected);break;
        }
    }
    std::ofstream(out/"COMPLETE")<<"Every pair validated; accepted DRAM/IO demand fills = 0.\n";
}
int main(int argc,char** argv)try {
    if(argc!=6)throw std::runtime_error("OUTPUT K near|mixed WINDOWS REPEATS");
    cpu_set_t set;CPU_ZERO(&set);CPU_SET(24,&set);
    if(sched_setaffinity(0,sizeof(set),&set))throw std::runtime_error("affinity");
    int k=std::stoi(argv[2]),repeats=std::stoi(argv[5]);size_t windows=std::stoull(argv[4]);
    std::string workload=argv[3];
    if(workload!="near"&&workload!="mixed")throw std::runtime_error("workload");
#define RUN(K) case K:run<K>(argv[1],workload,windows,repeats);break
    switch(k) {RUN(0);RUN(2);RUN(5);RUN(10);RUN(20);RUN(30);RUN(40);RUN(50);RUN(60);RUN(70);RUN(80);RUN(90);
        default:throw std::runtime_error("unsupported requested K");}
#undef RUN
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 2;}
