#include "backend/distributed/unstructured/fem/mars_mixed_overlap.hpp"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

using namespace mars::mixed;
static int checks=0;
static void require(bool ok,const char* reason) {
    ++checks;if(!ok) {std::fprintf(stderr,"FAIL: %s\n",reason);std::exit(1);}
}
static void near(double a,double b,double tolerance,const char* reason) {
    require(std::isfinite(a)&&std::isfinite(b)&&fabs(a-b)<=tolerance*std::max({1.,fabs(a),fabs(b)}),reason);
}
static Cell cell(int count,int begin) {
    Cell c;c.count=count;for(int k=0;k<count;++k)c.node[k]=begin+k;return c;
}
static void box(double* xyz,double x,double y,double z,double width=1) {
    for(int k=0;k<8;++k) {
        xyz[3*k]=x+width*.5*(1+sign(k,0));xyz[3*k+1]=y+width*.5*(1+sign(k,1));xyz[3*k+2]=z+.5*(1+sign(k,2));
    }
}
static double field(const double* x) {return 1+2*x[0]-3*x[1]+4*x[2];}
static std::vector<OverlapPatch> build(const Cell& a,const double* xa,int fa,const Cell& b,const double* xb,int fb) {
    std::vector<OverlapPatch> out;
    auto status=overlap_pair(a,xa,fa,b,xb,fb,{},[&](const OverlapPatch& p){out.push_back(p);});
    require(status==OverlapStatus::valid,"valid planar overlap");
    require(out.size()<=overlap_pair_capacity,"bounded per-pair patch count");
    bool topology[64]{};
    for(const auto& p:out) {
        require(p.topology_index>=0 && p.topology_index<64 && !topology[p.topology_index],"unique topological triangle-pair identity");
        topology[p.topology_index]=true;
        require(p.row>=0 && p.opposite_row>=0 && p.row!=p.opposite_row,"two distinct owner rows");
        require(dot(p.area,p.area)>0,"positive overlap area");
        for(int side=0;side<2;++side) {
            const Cell& c=side?b:a;const double* xyz=side?xb:xa;int offset=side?a.count:0;
            double sum=0,value=0,gradient[3]{};
            for(int k=0;k<c.count;++k) {
                sum+=p.shape[offset+k];value+=p.shape[offset+k]*field(xyz+3*k);
                for(int j=0;j<3;++j)gradient[j]+=p.gradient[3*(offset+k)+j]*field(xyz+3*k);
            }
            near(sum,1,1e-13,"constant interpolation");near(value,field(p.position),1e-12,"affine trace at shared physical point");
            near(gradient[0],2,1e-12,"affine gradient x");near(gradient[1],-3,1e-12,"affine gradient y");near(gradient[2],4,1e-12,"affine gradient z");
            require(p.inverse_height[side]>0,"positive penalty height");
        }
    }
    return out;
}
static double area(const std::vector<OverlapPatch>& patches) {
    double total=0;for(const auto& p:patches)total+=sqrt(dot(p.area,p.area));return total;
}
int main() {
    Cell a=cell(8,0),b=cell(8,8);double xa[24],xb[24];box(xa,0,0,0);box(xb,0,0,1);
    auto hh=build(a,xa,5,b,xb,4);near(area(hh),1,1e-13,"hex hex full overlap area");
    double row_area[16]{},moment=0;
    State state[16];Controls control;control.rho=1;control.mu=.1;
    for(int k=0;k<16;++k) {state[k].u[2]=.3;for(int j=0;j<3;++j)state[k].d[j]=1;}
    double residual[16]{};
    for(const auto& p:hh) {
        double weight=sqrt(dot(p.area,p.area));row_area[p.row]+=weight;row_area[p.opposite_row]+=weight;moment+=weight*field(p.position);
        auto f=interface_flux(p,state,control,0);residual[p.row]+=f.mass;residual[p.opposite_row]-=f.mass;
        near(f.mass,.3*weight,1e-14,"constant shared mass flux");
    }
    for(int k=0;k<4;++k) {near(row_area[4+k],.25,1e-13,"lower dual patch coverage");near(row_area[8+k],.25,1e-13,"upper dual patch coverage");}
    near(moment,4.5,1e-13,"affine integral at common centroids");
    double sum=0;for(double r:residual)sum+=r;near(sum,0,1e-15,"two-sided mass cancellation");
    auto reverse=build(b,xb,4,a,xa,5);near(area(reverse),area(hh),1e-13,"reversed pair same measure");
    for(const auto& p:reverse)require(p.area[2]<0,"reversed pair outward normal");

    double first[12]={0,0,1,1,0,1,0,1,1,0,0,2};
    double second[12]={1,0,1,1,1,1,0,1,1,1,0,2};
    Cell t0=cell(4,8),t1=cell(4,12);
    auto ht0=build(a,xa,5,t0,first,3),ht1=build(a,xa,5,t1,second,3);
    near(area(ht0)+area(ht1),1,1e-13,"two tet faces cover hex face");
    double lower[12]={0,0,1,0,1,1,1,0,1,0,0,0};
    auto tt=build(cell(4,0),lower,3,t0,first,3);near(area(tt),.5,1e-13,"tet tet area");

    double refined_area=0,coarse_rows[8]{};
    for(int i=0;i<2;++i)for(int j=0;j<2;++j) {
        box(xb,.5*i,.5*j,1,.5);auto patches=build(a,xa,5,b,xb,4);refined_area+=area(patches);
        for(const auto& p:patches)coarse_rows[p.row]+=sqrt(dot(p.area,p.area));
    }
    near(refined_area,1,1e-13,"nonmatching refinement covers face");
    for(int k=4;k<8;++k)near(coarse_rows[k],.25,1e-13,"refinement preserves coarse dual areas");
    box(xb,0,0,1);
    // Rotate and translate both cells to avoid relying on a coordinate plane.
    for(double* xyz:{xa,xb})for(int k=0;k<8;++k) {
        double x=xyz[3*k],y=xyz[3*k+1],z=xyz[3*k+2];
        xyz[3*k]=3+.6*x-.8*z;xyz[3*k+1]=-2+y;xyz[3*k+2]=4+.8*x+.6*z;
    }
    near(area(build(a,xa,5,b,xb,4)),1,1e-12,"rotated translated shared face measure");
    box(xa,0,0,0);
    box(xb,2,0,1);require(build(a,xa,5,b,xb,4).empty(),"disjoint coplanar faces emit no patches");
    box(xb,0,0,1);xb[2]+=.01;
    auto ignore=[](const OverlapPatch&){};
    require(overlap_pair(a,xa,5,b,xb,4,{},ignore)==OverlapStatus::invalid_geometry,"warped interface rejected");
    box(xb,0,0,1);
    require(overlap_pair(a,xa,5,b,xb,5,{},ignore)==OverlapStatus::invalid_geometry,"same normal rejected");
    require(overlap_pair(a,xa,5,a,xb,4,{},ignore)==OverlapStatus::invalid_input,"shared parent identities rejected");
    OverlapOptions invalid;invalid.distance_tolerance=-1;
    require(overlap_pair(a,xa,5,b,xb,4,invalid,ignore)==OverlapStatus::invalid_input,"invalid tolerance rejected");
    std::printf("PASS: %d shared mixed overlap checks\n",checks);
}
