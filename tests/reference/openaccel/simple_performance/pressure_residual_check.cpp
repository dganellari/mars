#include "../../../../backend/distributed/unstructured/fem/segregated/mars_segregated_pressure_capture.hpp"
#include "../../../../backend/distributed/unstructured/fem/segregated/mars_segregated_compensated_dot.hpp"
#include <iostream>
namespace frozen=mars::segregated::frozen;
using frozen::require;
static double down(double x) { return x==0?0:std::nextafter(x,0.); }
static double up(double x) { return std::nextafter(x,std::numeric_limits<double>::infinity()); }
struct NormInterval {
    double low=0,high=0;
    void add(double lower,double upper) {
        low=down(low+down(lower*lower)); high=up(high+up(upper*upper));
    }
    double lower() const { return down(std::sqrt(low)); }
    double upper() const { return up(std::sqrt(high)); }
};
int marker(const std::filesystem::path& directory,const char* schema) {
    std::ifstream in(directory/"complete"); std::string name,extra; int ranks=0;
    require(bool(in>>name>>ranks) && name==schema && ranks>0 && !(in>>extra)); return ranks;
}
int main(int argc,char** argv) {
    try {
        require(argc==3);
        const std::filesystem::path capture(argv[1]); const bool original=std::string(argv[2])=="-";
        const std::filesystem::path replay(argv[2]);
        const int ranks=marker(capture,"mars-pressure-capture-v1");
        if (!original) require(marker(replay,"mars-pressure-replay-v1")==ranks);
        const auto first=frozen::Part::read(capture/frozen::part_name(0));
        std::vector<double> owners(first.total),solution(first.total);
        std::uint64_t previous=0;
        for (int rank=0;rank<ranks;++rank) {
            const auto p=frozen::Part::read(capture/frozen::part_name(rank));
            require(p.rank==std::uint64_t(rank) && p.ranks==std::uint64_t(ranks) && p.first==previous && p.total==first.total
                && p.iteration==first.iteration && p.absolute==first.absolute && p.relative==first.relative && p.maximum==first.maximum);
            require(!p.solver_passed || !p.mars_passed);
            previous=p.last; std::vector<unsigned char> found(p.rows(),0);
            for (std::size_t local=0;local<p.nodes;++local) {
                const auto id=p.map[local];
                if (id>=std::int64_t(p.first) && id<std::int64_t(p.last)) {
                    require(!found[id-p.first]++ && std::isfinite(p.candidate[local])); owners[id]=p.candidate[local];
                }
            }
            require(std::all_of(found.begin(),found.end(),[](int value){return value==1;}));
            if (!original) {
                const auto file=replay/frozen::part_name(rank,".solution");
                require(std::filesystem::file_size(file)==p.rows()*sizeof(double));
                std::ifstream in(file,std::ios::binary);
                in.read(reinterpret_cast<char*>(solution.data()+p.first),std::streamsize(p.rows()*sizeof(double))); require(bool(in));
            }
        }
        require(previous==first.total);
        if (original) solution=owners;
        bool finite=std::all_of(solution.begin(),solution.end(),[](double x){return std::isfinite(x);});
        bool copies=true;
        NormInterval residual,rhs;
        for (int rank=0;rank<ranks;++rank) {
            const auto p=frozen::Part::read(capture/frozen::part_name(rank));
            for (std::size_t row=0;row<p.rows();++row) {
                mars::segregated::CompensatedDot dot; dot.product(p.rhs[row],1.);
                for (int k=p.offsets[row];k<p.offsets[row+1];++k) {
                    const auto local=p.columns[k]; const auto id=p.map[local];
                    copies=copies && p.candidate[local]==owners[id];
                    dot.product(-p.values[k],solution[id]);
                }
                const auto value=std::abs(dot.value()),error=dot.error_bound();
                finite=finite && std::isfinite(value) && std::isfinite(error);
                residual.add(down(std::max(0.,value-error)),up(value+error)); rhs.add(std::abs(p.rhs[row]),std::abs(p.rhs[row]));
            }
        }
        const double rel_low=down(first.relative*rhs.lower()),rel_high=up(first.relative*rhs.upper());
        const double limit_low=first.maximum?std::max(first.absolute,rel_low):down(first.absolute+rel_low);
        const double limit_high=first.maximum?std::max(first.absolute,rel_high):up(first.absolute+rel_high);
        finite=finite && std::isfinite(residual.upper()) && std::isfinite(rhs.upper()) && std::isfinite(limit_high);
        const bool passed=finite && residual.upper()<=limit_low;
        const bool failed=finite && residual.lower()>limit_high;
        std::cout<<std::boolalpha<<"{\"schema\":\"mars-pressure-residual-v1\",\"scope\":\"frozen_owned_rows_compensated_dot_with_evaluation_bound\","
            <<"\"capture_valid\":true,\"original_referenced_copies_equal_owners\":"<<copies
            <<",\"finite\":"<<finite<<",\"residual_passed\":"<<passed<<",\"residual_failed\":"<<failed
            <<",\"residual_inconclusive\":"<<(finite && !passed && !failed)<<"}\n";
        return 0;
    } catch (...) {
        std::cerr<<"ERROR: private pressure residual check failed\n"; return 1;
    }
}
