#pragma once
#include "backend/distributed/unstructured/fem/mars_mixed_overlap.hpp"

namespace mixedflow {
struct OverlapPair { int face[2]; };
struct OverlapIdentity { int row[2],pair,index; };
struct CountOverlap {
    int* count;
    MARS_FS_FN void operator()(const OverlapPatch&) const {++*count;}
};
struct StoreOverlap {
    Interface* patches;OverlapIdentity* identities;int pair,offset,capacity;int* count;
    MARS_FS_FN void operator()(const OverlapPatch& patch) const {
        const int index=(*count)++;
        if(index>=capacity) return;
        patches[offset+index]=patch;
        identities[offset+index]={{patch.row,patch.opposite_row},pair,patch.topology_index};
    }
};
struct BuildOverlap {
    View v;const OverlapPair* pairs;const int* offsets;int *counts,*valid;
    Interface* patches;OverlapIdentity* identities;
    MARS_FS_FN void operator()(int index) const {
        auto pair=pairs[index];const auto a=v.faces[pair.face[0]],b=v.faces[pair.face[1]];
        double xa[24],xb[24];v.coordinates(a.cell,xa);v.coordinates(b.cell,xb);int count=0;
        OverlapStatus status;
        if(!patches) status=overlap_pair(v.cells[a.cell],xa,a.ordinal,v.cells[b.cell],xb,b.ordinal,
                                      OverlapOptions{},CountOverlap{&count});
        else status=overlap_pair(v.cells[a.cell],xa,a.ordinal,v.cells[b.cell],xb,b.ordinal,
                                OverlapOptions{},StoreOverlap{patches,identities,index,offsets[index],counts[index],&count});
        valid[index]=status==OverlapStatus::valid && count<=overlap_pair_capacity && (!patches || count==counts[index]);
        if(!patches) counts[index]=count;
    }
};
struct CheckOverlapCoverage {
    View v;const int *offsets,*patches;int* valid;
    MARS_FS_FN void operator()(int q) const {
        const auto query=v.queries[q];const auto point=v.facepoint(query.face,query.sample,true);
        const double expected=sqrt(dot(point.area,point.area));double actual=0;
        for(int i=offsets[q];i<offsets[q+1];++i) {
            const auto& patch=v.interfaces[patches[i]];actual+=sqrt(dot(patch.area,patch.area));
        }
        valid[q]=std::isfinite(actual) && expected>0 && fabs(actual-expected)<=1e-8*expected;
    }
};

inline void Channel::build_conservative() {
    view_.conservative=true;view_.rank=rankId;
    const auto& input=input_;const auto& mesh=input.mesh;
    std::set<std::pair<int,int>> seen;std::vector<OverlapPair> pairs;
    for(auto query:input.queries) for(int i=query.start;i<query.start+query.count;++i) {
        int a=query.face,b=input.donors[i];
        if(input.faces[a].id>input.faces[b].id) std::swap(a,b);
        if(seen.insert({input.faces[a].id,input.faces[b].id}).second) pairs.push_back({{a,b}});
    }
    std::sort(pairs.begin(),pairs.end(),[&](auto a,auto b) {
        return std::make_pair(input.faces[a.face[0]].id,input.faces[a.face[1]].id)<
               std::make_pair(input.faces[b.face[0]].id,input.faces[b.face[1]].id);
    });
    Buffer<OverlapPair> d_pairs(pairs.begin(),pairs.end());Buffer<int> counts(pairs.size()),flags(pairs.size());
    each<device>(int(pairs.size()),BuildOverlap{view_,pointer(d_pairs),nullptr,pointer(counts),pointer(flags),nullptr,nullptr});
    ncflow::finish();ncflow::validFlags(flags,"mixed overlap geometry");
    // Counts and identities feed the existing host allocation/halo APIs; geometry
    // and the per-iteration flux remain on the device.
    const auto host_counts=host<int>(counts);std::vector<int> offsets{0};size_t total=0;
    for(int count:host_counts) {total+=size_t(count);if(total<=size_t(INT_MAX)) offsets.push_back(int(total));}
    check(total<=size_t(INT_MAX),"mixed overlap count range");
    Buffer<int> d_offsets(offsets.begin(),offsets.end());Buffer<OverlapIdentity> identities(total);
    interfaces_.resize(total);interface_mass_.resize(total);interface_fluxes_.resize(total);
    each<device>(int(pairs.size()),BuildOverlap{view_,pointer(d_pairs),pointer(d_offsets),pointer(counts),pointer(flags),pointer(interfaces_),pointer(identities)});
    ncflow::finish();ncflow::validFlags(flags,"mixed overlap construction");
    const auto metadata=host<OverlapIdentity>(identities);const auto owned=host<int>(owned_);
    std::vector<int> local_row(mesh.dofs.size(),-1);for(int r=0;r<n_;++r) local_row[owned[r]]=r;
    int local_max=-1,global_max=-1;for(auto face:input.faces) local_max=std::max(local_max,face.id);
    MPI_Allreduce(&local_max,&global_max,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    const uint64_t stride=uint64_t(global_max)+1,limit=std::numeric_limits<uint64_t>::max()/overlap_pair_capacity;
    check(stride>0 && stride<=limit/stride,"mixed overlap persistent identity range");
    std::vector<uint64_t> ids;std::vector<int> owners;
    std::vector<std::vector<int>> by_row(n_),signs(n_);
    std::map<std::pair<int,int>,std::vector<int>> coverage;
    for(size_t q=0;q<metadata.size();++q) {
        const auto meta=metadata[q];const auto pair=pairs[meta.pair];
        ids.push_back((uint64_t(input.faces[pair.face[0]].id)*stride+input.faces[pair.face[1]].id)*overlap_pair_capacity+meta.index);
        owners.push_back(owner(mesh.dofs[meta.row[0]],mesh.globalNodes));
        for(int side=0;side<2;++side) {
            int r=local_row[meta.row[side]];if(r<0) continue;
            by_row[r].push_back(int(q));signs[r].push_back(side?-1:1);
            coverage[{pair.face[side],meta.row[side]}].push_back(int(q));
        }
    }
    ncflow::incidence(by_row,qo_,qi_);std::vector<int> flat_signs;
    for(const auto& list:signs) flat_signs.insert(flat_signs.end(),list.begin(),list.end());
    qs_=Buffer<int>(flat_signs.begin(),flat_signs.end());patch_owner_=Buffer<int>(owners.begin(),owners.end());
    view_.interfaces=pointer(interfaces_);view_.interface_mass=pointer(interface_mass_);view_.interface_fluxes=pointer(interface_fluxes_);
    view_.qo=pointer(qo_);view_.qi=pointer(qi_);view_.qs=pointer(qs_);view_.patch_owner=pointer(patch_owner_);
    std::vector<std::vector<int>> query_patches;for(auto query:input.queries) query_patches.push_back(coverage[{query.face,query.row}]);
    Buffer<int> coverage_offsets,coverage_ids;ncflow::incidence(query_patches,coverage_offsets,coverage_ids);
    flags.resize(input.queries.size());each<device>(int(input.queries.size()),CheckOverlapCoverage{view_,pointer(coverage_offsets),pointer(coverage_ids),pointer(flags)});
    ncflow::finish();ncflow::validFlags(flags,"mixed overlap surface coverage");
    check(flux_halo_.buildFromHost(ids,owners,mesh.candidatePeers),"mixed overlap flux ownership");
    check(mass_halo_.buildFromHost(ids,owners,mesh.candidatePeers),"mixed overlap mass ownership");
}
}
