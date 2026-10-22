#pragma once
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <cstring>
#include <x86intrin.h>
struct Event { const char*name;uint32_t type;uint64_t config; };
inline std::array<Event,2> events(const std::string& name){
    if(name=="l1d")return {{{"l1d_access",PERF_TYPE_HW_CACHE,PERF_COUNT_HW_CACHE_L1D},
        {"l1d_miss",PERF_TYPE_HW_CACHE,PERF_COUNT_HW_CACHE_L1D|(uint64_t(PERF_COUNT_HW_CACHE_RESULT_MISS)<<16)}}};
    if(name=="l2")return {{{"amd_l2_requests",PERF_TYPE_HARDWARE,PERF_COUNT_HW_CACHE_REFERENCES},
        {"amd_l2_ic_dc_miss",PERF_TYPE_HARDWARE,PERF_COUNT_HW_CACHE_MISSES}}};
    if(name=="dram")return {{{"zen4_demand_fills_all",PERF_TYPE_RAW,0xff43},
        {"zen4_demand_fills_dram_io",PERF_TYPE_RAW,0x4843}}};
    if(name=="frontend")return {{{"l1i_miss",PERF_TYPE_HW_CACHE,PERF_COUNT_HW_CACHE_L1I|(uint64_t(PERF_COUNT_HW_CACHE_RESULT_MISS)<<16)},
        {"dtlb_walk",PERF_TYPE_HW_CACHE,PERF_COUNT_HW_CACHE_DTLB|(uint64_t(PERF_COUNT_HW_CACHE_RESULT_MISS)<<16)}}};
    if(name=="work")return {{{"user_cycles",PERF_TYPE_HARDWARE,PERF_COUNT_HW_CPU_CYCLES},{"user_instructions",PERF_TYPE_HARDWARE,PERF_COUNT_HW_INSTRUCTIONS}}};
    throw std::runtime_error("event group");
}
struct PerfPair{
    int fd[2]{-1,-1};perf_event_mmap_page* page[2]{};std::array<Event,2> spec;
    uint64_t enabled=0,running=0,total[2]{};
    explicit PerfPair(const std::string& group):spec(events(group)){
        for(int i=0;i<2;++i){perf_event_attr a{};a.size=sizeof(a);a.type=spec[i].type;a.config=spec[i].config;
            a.disabled=i==0;a.pinned=i==0;a.exclude_kernel=1;a.exclude_hv=1;
            a.read_format=PERF_FORMAT_GROUP|PERF_FORMAT_TOTAL_TIME_ENABLED|PERF_FORMAT_TOTAL_TIME_RUNNING;
            fd[i]=syscall(SYS_perf_event_open,&a,0,-1,i?fd[0]:-1,0);
            if(fd[i]<0)throw std::runtime_error(std::string(spec[i].name)+": "+std::strerror(errno));
            page[i]=(perf_event_mmap_page*)mmap(nullptr,4096,PROT_READ,MAP_SHARED,fd[i],0);
            if(page[i]==MAP_FAILED||!page[i]->cap_user_rdpmc)throw std::runtime_error("user RDPMC unavailable");
        }
    }
    void start(){if(ioctl(fd[0],PERF_EVENT_IOC_RESET,PERF_IOC_FLAG_GROUP)||ioctl(fd[0],PERF_EVENT_IOC_ENABLE,PERF_IOC_FLAG_GROUP))throw std::runtime_error("perf start");}
    __attribute__((always_inline)) static inline uint64_t fast(perf_event_mmap_page* p){
        uint32_t seq,index,width;int64_t offset,pmc;
        for(;;){
            seq=__atomic_load_n(&p->lock,__ATOMIC_ACQUIRE);asm volatile("" ::: "memory");
            if(seq&1)continue;
            index=p->index;offset=p->offset;width=p->pmc_width;
            if(!index)throw std::runtime_error("pinned event not scheduled");
            pmc=(int64_t)(__rdpmc(index-1)<<(64-width))>>(64-width);offset+=pmc;
            asm volatile("" ::: "memory");
            if(__atomic_load_n(&p->lock,__ATOMIC_RELAXED)==seq)break;
        }
        return uint64_t(offset);
    }
    __attribute__((always_inline)) inline std::array<uint64_t,2> snapshot(){
        _mm_lfence();auto a=fast(page[0]),b=fast(page[1]);_mm_lfence();return {a,b};
    }
    void stop(){
        if(ioctl(fd[0],PERF_EVENT_IOC_DISABLE,PERF_IOC_FLAG_GROUP))throw std::runtime_error("perf stop");
        uint64_t v[5]{};if(read(fd[0],v,sizeof(v))!=sizeof(v)||v[0]!=2)throw std::runtime_error("perf group read");
        enabled=v[1];running=v[2];total[0]=v[3];total[1]=v[4];
        if(!enabled||double(running)/enabled<.99)throw std::runtime_error("perf multiplexed / unscheduled");
    }
    ~PerfPair(){for(int i=0;i<2;++i){if(page[i]&&page[i]!=MAP_FAILED)munmap(page[i],4096);if(fd[i]>=0)close(fd[i]);}}
};
