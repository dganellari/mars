#include "mars_segregated_assembly.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>
using namespace mars::segregated;
int checks=0;
void require(bool value) { ++checks; if (!value) throw std::runtime_error("assembly algebra check failed"); }
void close(double a,double b) { require(std::isfinite(a) && std::abs(a-b)<=1e-13*std::max(1.,std::abs(b))); }

template<int Components> void shared_node_fixture() {
    int offsets[]{0,2,5,7},columns[]{0,1,0,1,2,1,2};
    std::vector<double> values(7*Components*Components,0),rhs(3*Components,0);
    std::vector<double> dense(9*Components*Components,0),expected_rhs(3*Components,0);
    BlockCsrView<Components> matrix{3,offsets,columns,values.data(),rhs.data()};
    for (int cell=0;cell<2;++cell) {
        // Reverse the second local node order to exercise both mappings at the shared node.
        int nodes[]{cell==0?0:2,1};
        std::array<double,4*Components*Components> lhs;
        std::array<double,2*Components> local_rhs;
        for (int r=0;r<2*Components;++r) {
            local_rhs[r]=(cell+1)*(.5+r);
            expected_rhs[nodes[r/Components]*Components+r%Components]+=local_rhs[r];
            for (int c=0;c<2*Components;++c) {
                lhs[r*2*Components+c]=1+cell+2*r+.25*c;
                const int row=nodes[r/Components]*Components+r%Components;
                const int col=nodes[c/Components]*Components+c%Components;
                dense[row*3*Components+col]+=lhs[r*2*Components+c];
            }
        }
        require(scatter_block(matrix,nodes,2,lhs.data(),local_rhs.data()));
    }
    for (int row=0;row<3;++row) {
        for (int i=0;i<Components;++i) close(rhs[row*Components+i],expected_rhs[row*Components+i]);
        for (int k=offsets[row];k<offsets[row+1];++k)
            for (int i=0;i<Components;++i) for (int j=0;j<Components;++j)
                close(values[k*Components*Components+i*Components+j],
                      dense[(row*Components+i)*3*Components+columns[k]*Components+j]);
    }
    auto old_values=values,old_rhs=rhs;
    double local_lhs[16*Components*Components]{},local_rhs[4*Components]{};
    const int missing[]{0,2},invalid[]{-1},outside[]{3};
    require(!scatter_block(matrix,missing,2,local_lhs,local_rhs));
    require(!scatter_block(matrix,invalid,1,local_lhs,local_rhs));
    require(!scatter_block(matrix,outside,1,local_lhs,local_rhs));
    require(!scatter_block(matrix,missing,0,local_lhs,local_rhs));
    require(values==old_values && rhs==old_rhs); // No partial write before an absent zero entry is rejected.
    if constexpr (Components==3) {
        double d[3],dt[3]; const int diagonal=3;
        require(finish_momentum_row(matrix,1,2.,.5,.75,false,d,dt));
        for (int k=0;k<7;++k) for (int i=0;i<3;++i) for (int j=0;j<3;++j)
            close(values[9*k+3*i+j],old_values[9*k+3*i+j]*(k==diagonal && i==j?2.:1.));
        for (int i=0;i<3;++i) {
            close(d[i],2./(2*old_values[9*diagonal+3*i+i]+std::numeric_limits<double>::epsilon()));
            close(rhs[3+i],.75*old_rhs[3+i]); close(dt[i],0);
        }
        require(!finish_momentum_row(matrix,1,2.,0.,.75,false,d,dt));
        require(!finish_momentum_row(matrix,3,2.,.5,.75,false,d,dt));
    }
}
int main() {
    try {
        shared_node_fixture<1>(); shared_node_fixture<3>();
        // Opposite-node outlet derivative: row is the face node, column the opposite node.
        int offsets[]{0,2,4},columns[]{0,1,0,1},nodes[]{0,1};
        double values[4]{},rhs[2]{},local_lhs[]{0,3,0,0},local_rhs[]{-2,2};
        BlockCsrView<1> matrix{2,offsets,columns,values,rhs};
        require(scatter_block(matrix,nodes,2,local_lhs,local_rhs));
        close(values[1],3); close(values[2],0); close(rhs[0]+rhs[1],0);
        close(values[0]+values[1],3); // The compact outlet term anchors constants without a pinned row.
        std::cout<<"PASS: "<<checks<<" independent block assembly checks\n";
    } catch (const std::exception& e) { std::cerr<<"FAIL: "<<e.what()<<'\n'; return 1; }
}
