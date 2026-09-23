#pragma once
#include "export_writer.hpp"
#include <mpi.h>
#include <cstdlib>
#include <iostream>

namespace mars_reference {
struct UpdateChunk {
    const double* data = nullptr;
    std::size_t size = 1;
    double scalar = 0;
    UpdateChunk(double value) : scalar(value) {}
    UpdateChunk(const double* values, std::size_t n) : data(values), size(n) {}
};
class UpdateSession;
inline UpdateSession* current_update = nullptr;
class UpdateSession {
    MPI_Comm comm_;
    std::ofstream out_;
    int phase_ = -1;
    std::size_t count_ = 0;
    template<class F> void checked(F&& action) {
        try { action(); }
        catch(const std::exception& e) {
            std::cerr << "Update export failed: " << e.what() << '\n';
            MPI_Abort(comm_,1); std::abort();
        }
    }
    void array(std::initializer_list<UpdateChunk> chunks) {
        out_ << '['; bool first=true;
        for(const auto& chunk:chunks) for(std::size_t i=0;i<chunk.size;++i) {
            const double value=chunk.data ? chunk.data[i] : chunk.scalar;
            require(std::isfinite(value),"nonfinite update value");
            if(!first) out_<<',';
            first=false; out_<<value;
        }
        out_<<']';
    }
public:
    explicit UpdateSession(MPI_Comm comm, bool supported = true) : comm_(comm) {
        const char* root=std::getenv("MARS_OPENACCEL_EXPORT_DIR");
        if(!root || !*root) return;
        checked([&] {
            require(supported,"updates require steady 3D flow with one pressure subiteration");
            const char* fixture=std::getenv("MARS_OPENACCEL_PUBLIC_FIXTURE");
            require(fixture && std::string(fixture)=="public_channel","updates require the public channel");
            require(!current_update,"nested update export");
            int ranks=0;
            require(MPI_Comm_size(comm,&ranks)==MPI_SUCCESS && ranks==1,"update capture requires one rank");
            static int iteration=0;
            const auto directory=std::filesystem::path(root)/"updates";
            require(std::filesystem::is_directory(root),"create fresh export root first");
            std::filesystem::create_directory(directory);
            const auto path=directory/("iteration"+std::to_string(++iteration)+".jsonl");
            require(!std::filesystem::exists(path),"update export exists");
            out_.exceptions(std::ios::failbit|std::ios::badbit); out_.open(path);
            out_.imbue(std::locale::classic());
            out_<<std::setprecision(17)<<"{\"kind\":\"header\",\"schema\":1,\"fixture\":\"public_channel\","
                <<"\"producer\":\"openaccel\",\"reference_revision\":\""<<reference_revision
                <<"\",\"solver_revision\":\""<<solver_revision<<"\",\"iteration\":"<<iteration<<",\"ranks\":1}\n";
            current_update=this;
        });
    }
    void phase(int value) {
        if(!out_.is_open()) return;
        checked([&] {
            require(value==phase_+1,"unexpected SIMPLE update order or repeated pressure subiteration");
            phase_=value; out_<<"{\"kind\":\"phase\",\"value\":"<<value<<"}\n";
        });
    }
    void record(int stage, std::uint64_t entity, int sample,
                std::initializer_list<UpdateChunk> inputs, std::initializer_list<UpdateChunk> outputs) {
        checked([&] {
            require(entity>0 && sample>=0,"invalid update identity");
            out_<<"{\"kind\":\"update\",\"stage\":"<<stage<<",\"entity\":"<<entity<<",\"sample\":"<<sample
                <<",\"phase\":"<<phase_<<",\"inputs\":";
            array(inputs); out_<<",\"outputs\":"; array(outputs); out_<<"}\n"; ++count_;
        });
    }
    int phase() const { return phase_; }
    ~UpdateSession() {
        if(out_.is_open()) checked([&] {
            require(phase_==7,"incomplete SIMPLE update sequence");
            out_<<"{\"kind\":\"end\",\"records\":"<<count_<<"}\n";out_.close();current_update=nullptr;
        });
    }
};
inline bool update_active(int phase=-1) {
    return current_update && (phase<0 || current_update->phase()==phase);
}
inline void update_phase(int phase) { if(current_update) current_update->phase(phase); }
inline void update_record(int stage,std::uint64_t entity,int sample,
                          std::initializer_list<UpdateChunk> inputs,std::initializer_list<UpdateChunk> outputs) {
    if(current_update) current_update->record(stage,entity,sample,inputs,outputs);
}
} // namespace mars_reference
