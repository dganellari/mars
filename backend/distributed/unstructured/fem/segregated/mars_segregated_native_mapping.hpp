#pragma once
#include "mars_segregated_simple.hpp"
#ifdef __CUDACC__
#define MARS_NATIVE_HD __host__ __device__
#else
#define MARS_NATIVE_HD
#endif

namespace mars::segregated {
// Sorted node identities survive element reordering and local vertex permutations.
struct SimpleCellKey {
    int nodes[4]{};
    MARS_NATIVE_HD bool operator<(const SimpleCellKey& other) const {
        for (int j=0;j<4;++j) {
            if (nodes[j]<other.nodes[j]) return true;
            if (nodes[j]>other.nodes[j]) return false;
        }
        return false;
    }
    MARS_NATIVE_HD bool operator==(const SimpleCellKey& other) const {
        for (int j=0;j<4;++j) if (nodes[j]!=other.nodes[j]) return false;
        return true;
    }
};
MARS_NATIVE_HD inline SimpleCellKey simple_cell_key(const int* nodes) {
    SimpleCellKey key;
    for (int j=0;j<4;++j) key.nodes[j]=nodes[j];
    for (int i=1;i<4;++i) for (int j=i;j>0 && key.nodes[j]<key.nodes[j-1];--j) {
        int t=key.nodes[j]; key.nodes[j]=key.nodes[j-1]; key.nodes[j-1]=t;
    }
    return key;
}
MARS_NATIVE_HD inline int simple_native_face(const int* source_face,const int* cell) {
    int match=-1;
    for (int f=0;f<4;++f) {
        int count=0;
        for (int j=0;j<3;++j) for (int k=0;k<3;++k)
            count+=source_face[j]==cell[tet_face_node(f,k)];
        if (count==3) { if (match>=0) return -1; match=f; }
    }
    return match;
}
template<class KeyType>
MARS_NATIVE_HD int simple_find_key(const KeyType* keys,int size,const KeyType& key) {
    int lo=0,hi=size;
    while (lo<hi) { int mid=lo+(hi-lo)/2; if (keys[mid]<key) lo=mid+1; else hi=mid; }
    return lo<size && keys[lo]==key?lo:-1;
}
}
#undef MARS_NATIVE_HD
