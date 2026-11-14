#include "examples/distributed/unstructured/mars_mixed_nonconformal_runtime.hpp"
#include <array>
#include <cstdio>
#include <stdexcept>
using namespace mixedflow;
static int checks;
static void require(bool condition,const char* message) {
    ++checks;if(!condition) throw std::runtime_error(message);
}
static void near(double a,double b,const char* message) {require(std::abs(a-b)<1e-13,message);}
static void face_scores() {
    for(int count:{4,8}) {
        double xyz[24]{};
        if(count==4) for(int k=1;k<4;++k) xyz[3*k+k-1]=1;
        else for(int k=0;k<8;++k) for(int j=0;j<3;++j) xyz[3*k+j]=sign(k,j);
        for(int face=0;face<face_count(count);++face) for(double r:{-.8,-.2,0.,.6,1.}) {
            double reference[3]{};
            if(count==4) {
                double shape[3]={.2,.4+.1*r,.4-.1*r};
                for(int k=0;k<3;++k) {int n=face_node(count,face,k);if(n) reference[n-1]=shape[k];}
                Point p;require(evaluate(count,xyz,reference,p),"triangle parent evaluation");
                // Nalu uses max(-3*(lambda0-1/3),-3*(lambda1-1/3),
                // 3*(lambda0+lambda1-2/3)) at zero normal distance.
                double x=shape[0]-1./3,y=shape[1]-1./3;
                near(donor_face_score(count,face,p),std::max({-3*x,-3*y,3*(x+y)}),"triangle Nalu distance");
            } else {
                const double s=-.37;
                const double shape[4]={.25*(1-r)*(1-s),.25*(1+r)*(1-s),.25*(1+r)*(1+s),.25*(1-r)*(1+s)};
                for(int k=0;k<4;++k) for(int j=0;j<3;++j) reference[j]+=shape[k]*sign(face_node(count,face,k),j);
                Point p;require(evaluate(count,xyz,reference,p),"quad parent evaluation");
                near(donor_face_score(count,face,p),std::max(std::abs(r),std::abs(s)),"quad Nalu distance");
            }
        }
    }
}
static void candidate_selection() {
    // Current face is below z=0. Two valid donors cover its first vertex:
    // one at a donor vertex, one at the donor centroid.
    double xyz[]={0,0,0, 0,1,0, 1,0,0, 0,0,-1,
                        0,0,0, 2,0,0, 0,2,0, 0,0,1,
                        -1,-1,0, 1,-1,0, 0,2,0, 0,0,1};
    Cell cells[3]{};Face faces[3]{};
    for(int e=0;e<3;++e) {
        cells[e].count=4;for(int k=0;k<4;++k) cells[e].node[k]=4*e+k;
        int ordinal=-1;
        for(int f=0;f<4;++f) {
            bool plane=true;for(int k=0;k<3;++k) plane=plane && xyz[3*(4*e+face_node(4,f,k))+2]==0;
            if(plane) ordinal=f;
        }
        require(ordinal>=0,"interface face found");faces[e]={e,ordinal,e==0?3:4,e};
    }
    int sample=-1;for(int k=0;k<3;++k) if(face_node(4,faces[0].ordinal,k)==0) sample=k;
    require(sample>=0,"current vertex found");
    Query query{0,sample,0,2,0};int candidates[2]={1,2},valid[1];Interface result[1];
    View v{};v.cells=cells;v.faces=faces;v.queries=&query;v.donors=candidates;v.xyz=xyz;v.interfaces=result;
    for(int repeat=0;repeat<2;++repeat) {
        BuildInterface{v,valid,DonorSelection::persistent_id}(0);
        require(valid[0] && result[0].node[4]==4,"legacy metric keeps lowest-ID donor");
        BuildInterface{v,valid,DonorSelection::openaccel_face_score}(0);
        require(valid[0] && result[0].node[4]==8,"face score selects donor interior over donor edge");
        std::swap(candidates[0],candidates[1]);
    }
    // Equal geometry selects persistent identity regardless of traversal order.
    for(int k=0;k<12;++k) xyz[12+k]=xyz[24+k];
    for(int repeat=0;repeat<2;++repeat) {
        BuildInterface{v,valid,DonorSelection::openaccel_face_score}(0);
        require(valid[0] && result[0].node[4]==4,"equal score selects persistent donor identity");
        std::swap(candidates[0],candidates[1]);
    }
}
static void outlet_storage() {
    Cell cell{};cell.count=8;double xyz[24]{};State state[8]{};
    for(int k=0;k<8;++k) {cell.node[k]=k;for(int j=0;j<3;++j) xyz[3*k+j]=.5*(1+sign(k,j));state[k].u[0]=3;}
    Face face{0,1,1,0};state[face_node(8,1,0)].u[0]=-5;
    double mass[4]{},expected[4]{},net=0;int closed=0;View v{};
    v.cells=&cell;v.faces=&face;v.xyz=xyz;v.state=state;v.face_mass=mass;v.closed=&closed;v.completed=1;
    bool negative=false;
    for(int q=0;q<4;++q) {expected[q]=boundary_flux(cell,1,v.facepoint(0,q),state,xyz,v.controls,1,0,false,true).mass;net+=expected[q];negative|=expected[q]<0;}
    require(net>0 && negative,"open outlet fixture has mixed-sign point fluxes");
    UpdateBoundaryMass{v}(0);require(!closed,"positive net outlet remains open");
    for(int q=0;q<4;++q) near(mass[q],expected[q],"stored outlet flux equals pressure equation flux");
    for(auto& s:state) s.u[0]=-1;
    UpdateBoundaryMass{v}(0);require(closed,"net backflow closes entire face");
    for(double value:mass) near(value,0,"closed outlet has zero mass");
}
static void missing_column() {
    int64_t ids[3]={0,1,2},columns[2]={0,1};int offsets[2]={0,2},slots[2]={0,0},valid=1;
    double values[2]{};View v{};v.ids=ids;v.coefficient_slots=slots;
    InterfaceCsrView<double> matrix{};matrix.columns=columns;matrix.offsets=offsets;matrix.values=values;
    Assemble assembly{v,matrix,nullptr,&valid,false};
    assembly.add(0,1,0,0,2);near(values[0],2,"declared omitted column folds to diagonal");
    assembly.add(0,2,0,0,3);require(!valid,"folding must not hide a genuinely missing column");
    near(values[0],2,"missing coefficient is not silently folded");
}
int main() {
    face_scores();candidate_selection();outlet_storage();missing_column();
    std::printf("PASS: %d coincident-face score and production donor checks\n",checks);
}
