#pragma once
#include "perf_pair.hpp"
struct PerfGate{
    int fd[3]{-1,-1,-1};perf_event_mmap_page* page[3]{};
    uint64_t enabled=0,running=0,total[3]{};
    PerfGate(){
        std::array<Event,3> spec{{{"user_cycles",PERF_TYPE_HARDWARE,PERF_COUNT_HW_CPU_CYCLES},
            {"amd_l2_ic_dc_miss",PERF_TYPE_HARDWARE,PERF_COUNT_HW_CACHE_MISSES},
            {"zen4_demand_fills_dram_io",PERF_TYPE_RAW,0x4843}}};
        for(int i=0;i<3;++i){perf_event_attr a{};a.size=sizeof(a);a.type=spec[i].type;a.config=spec[i].config;
            a.disabled=i==0;a.pinned=i==0;a.exclude_kernel=1;a.exclude_hv=1;
            a.read_format=PERF_FORMAT_GROUP|PERF_FORMAT_TOTAL_TIME_ENABLED|PERF_FORMAT_TOTAL_TIME_RUNNING;
            fd[i]=syscall(SYS_perf_event_open,&a,0,-1,i?fd[0]:-1,0);
            if(fd[i]<0)throw std::runtime_error(std::string(spec[i].name)+": "+std::strerror(errno));
            page[i]=(perf_event_mmap_page*)mmap(nullptr,4096,PROT_READ,MAP_SHARED,fd[i],0);
            if(page[i]==MAP_FAILED||!page[i]->cap_user_rdpmc)throw std::runtime_error("user RDPMC unavailable");
        }
    }
    void start(){if(ioctl(fd[0],PERF_EVENT_IOC_RESET,PERF_IOC_FLAG_GROUP)||ioctl(fd[0],PERF_EVENT_IOC_ENABLE,PERF_IOC_FLAG_GROUP))throw std::runtime_error("gate start");}
    __attribute__((always_inline)) inline std::array<uint64_t,3> snapshot(){
        _mm_lfence();auto a=PerfPair::fast(page[0]),b=PerfPair::fast(page[1]),c=PerfPair::fast(page[2]);_mm_lfence();return {a,b,c};
    }
    void stop(){
        if(ioctl(fd[0],PERF_EVENT_IOC_DISABLE,PERF_IOC_FLAG_GROUP))throw std::runtime_error("gate stop");
        uint64_t v[6]{};if(read(fd[0],v,sizeof(v))!=sizeof(v)||v[0]!=3)throw std::runtime_error("gate group read");
        enabled=v[1];running=v[2];for(int i=0;i<3;++i)total[i]=v[3+i];
        if(!enabled||double(running)/enabled<.99)throw std::runtime_error("gate multiplexed / unscheduled");
    }
    ~PerfGate(){for(int i=0;i<3;++i){if(page[i]&&page[i]!=MAP_FAILED)munmap(page[i],4096);if(fd[i]>=0)close(fd[i]);}}
};
