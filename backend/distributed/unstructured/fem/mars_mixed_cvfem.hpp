#pragma once
#include "segregated/mars_segregated_geometry.hpp"

#if defined(__CUDACC__) || defined(__HIPCC__)
#define MARS_MIXED_HD __host__ __device__
#else
#define MARS_MIXED_HD
#endif

namespace mars::mixed {
struct Cell { int node[8]{}, count{}; };
struct Point {
    double shape[8]{}, gradient[24]{}, area[3]{}, position[3]{};
    int left{}, right{};
};
struct Geometry { Point interior[12]{}; double volume[8]{}; int samples{}; };
MARS_MIXED_HD inline int face_size(int count) { return count==8 ? 4 : 3; }
MARS_MIXED_HD inline int face_count(int count) { return count==8 ? 6 : 4; }
MARS_MIXED_HD inline int face_node(int count,int face,int k) {
    constexpr int hex[24]={0,1,5,4, 1,2,6,5, 2,3,7,6, 0,4,7,3, 0,3,2,1, 4,5,6,7};
    return count==8 ? hex[4*face+k] : segregated::tet_face_node(face,k);
}
MARS_MIXED_HD inline int sign(int node,int axis) {
    constexpr int s[24]={-1,-1,-1, 1,-1,-1, 1,1,-1, -1,1,-1,
                        -1,-1,1, 1,-1,1, 1,1,1, -1,1,1};
    return s[3*node+axis];
}
MARS_MIXED_HD inline void cross(const double* a,const double* b,double* c) {
    for(int j=0;j<3;++j) c[j]=a[(j+1)%3]*b[(j+2)%3]-a[(j+2)%3]*b[(j+1)%3];
}
MARS_MIXED_HD inline double dot(const double* a,const double* b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
// Reference hex coordinates are [-1,1]^3. The tetrahedron uses barycentric r,s,t.
MARS_MIXED_HD inline bool evaluate(int count,const double* xyz,const double* reference,Point& p,double* jacobian=nullptr,double* determinant=nullptr) {
    double derivative[24]{},j[9]{},cofactor[9];
    for(int k=0;k<8;++k) p.shape[k]=0;
    if(count==8) {
        for(int k=0;k<8;++k) {
            p.shape[k]=.125;
            for(int a=0;a<3;++a) p.shape[k]*=1+sign(k,a)*reference[a];
            for(int a=0;a<3;++a) {
                derivative[3*k+a]=.125*sign(k,a);
                for(int b=0;b<3;++b) if(b!=a) derivative[3*k+a]*=1+sign(k,b)*reference[b];
            }
        }
    } else if(count==4) {
        p.shape[0]=1-reference[0]-reference[1]-reference[2];
        for(int a=0;a<3;++a) {p.shape[a+1]=reference[a];derivative[a]=-1;derivative[3*(a+1)+a]=1;}
    } else return false;
    for(int a=0;a<3;++a) {
        p.position[a]=0;
        for(int k=0;k<count;++k) {
            p.position[a]+=p.shape[k]*xyz[3*k+a];
            for(int b=0;b<3;++b) j[3*b+a]+=xyz[3*k+a]*derivative[3*k+b];
        }
    }
    cross(j+3,j+6,cofactor);cross(j+6,j,cofactor+3);cross(j,j+3,cofactor+6);
    double det=dot(j,cofactor);
    if(!(det>0) || !segregated::geometry_finite(det)) return false;
    for(int k=0;k<count;++k) for(int a=0;a<3;++a) {
        p.gradient[3*k+a]=0;
        for(int b=0;b<3;++b) p.gradient[3*k+a]+=derivative[3*k+b]*cofactor[3*b+a]/det;
    }
    if(jacobian) for(int a=0;a<9;++a) jacobian[a]=j[a];
    if(determinant) *determinant=det;
    return true;
}
MARS_MIXED_HD inline bool geometry(int count,const double* xyz,Geometry& g,bool shifted_gradients=false) {
    g={};
    if(count==4) {
        segregated::TetGeometry<double> tet;
        if(!segregated::tet_geometry(xyz,tet)) return false;
        g.samples=6;
        for(int k=0;k<4;++k) g.volume[k]=tet.volume/4;
        for(int q=0;q<6;++q) {
            auto& p=g.interior[q]; double shape[4];segregated::tet_sample_shape(q,false,shape);
            if(!evaluate(4,xyz,shape+1,p)) return false;
            p.left=segregated::tet_edge_node(q,0);p.right=segregated::tet_edge_node(q,1);
            for(int a=0;a<3;++a) p.area[a]=tet.area[3*q+a];
        }
    } else if(count==8) {
        g.samples=12;int q=0;
        for(int axis=0;axis<3;++axis) for(int k=0;k<8;++k) if(sign(k,axis)<0) {
            auto& p=g.interior[q++];double r[3],j[9];
            for(int a=0;a<3;++a) r[a]=a==axis?0:.5*sign(k,a);
            if(!evaluate(8,xyz,r,p,j)) return false;
            p.left=k;
            for(int l=0;l<8;++l) {
                bool match=sign(l,axis)>0;
                for(int a=0;a<3;++a) if(a!=axis) match=match && sign(l,a)==sign(k,a);
                if(match) p.right=l;
            }
            cross(j+3*((axis+1)%3),j+3*((axis+2)%3),p.area);
            if(shifted_gradients) {
                // Flow derivatives use the edge midpoint; interpolation and areas stay at the subface point.
                for(int a=0;a<3;++a) r[a]=.5*(sign(p.left,a)+sign(p.right,a));
                Point derivative;if(!evaluate(8,xyz,r,derivative)) return false;
                for(int k=0;k<24;++k) p.gradient[k]=derivative.gradient[k];
            }
        }
        // Integrate each median-dual subvolume, including trilinear Jacobian variation.
        constexpr double abscissa=.57735026918962576451;
        for(int k=0;k<8;++k) for(int q=0;q<8;++q) {
            double r[3],det;Point p;
            for(int a=0;a<3;++a) r[a]=.5*sign(k,a)+.5*abscissa*sign(q,a);
            if(!evaluate(8,xyz,r,p,nullptr,&det)) return false;
            g.volume[k]+=det*.125;
        }
    } else return false;
    return true;
}
MARS_MIXED_HD inline bool face_point(int count,const double* xyz,int face,int sample,bool shifted,Point& p,bool shifted_gradients=false) {
    if(face<0 || face>=face_count(count) || sample<0 || sample>=face_size(count)) return false;
    int node=face_node(count,face,sample);
    if(count==4) {
        double shape[3],r[3]{};segregated::tri_sample_shape(sample,shifted,shape);
        for(int k=0;k<3;++k) {int c=face_node(count,face,k);if(c) r[c-1]=shape[k];}
        if(!evaluate(count,xyz,r,p)) return false;
        segregated::TetGeometry<double> g;
        if(!segregated::tet_geometry(xyz,g)) return false;
        segregated::tet_boundary_area(g,face,p.area);
    } else {
        int axis=-1,normal=0;
        for(int a=0;a<3;++a) {
            bool same=true;for(int k=1;k<4;++k) same=same && sign(face_node(8,face,k),a)==sign(face_node(8,face,0),a);
            if(same) {axis=a;normal=sign(node,a);}
        }
        if(axis<0) return false;
        double r[3],j[9];
        for(int a=0;a<3;++a) r[a]=a==axis?normal:(shifted?1.:.5)*sign(node,a);
        if(!evaluate(count,xyz,r,p)) return false;
        // Areas use the ordinary subface point even when interpolation is shifted.
        for(int a=0;a<3;++a) r[a]=a==axis?normal:.5*sign(node,a);
        Point area_point;if(!evaluate(count,xyz,r,area_point,j)) return false;
        cross(j+3*((axis+1)%3),j+3*((axis+2)%3),p.area);
        for(int a=0;a<3;++a) p.area[a]*=normal;
        if(shifted_gradients && !shifted) {
            // The boundary derivative is evaluated at the nearest vertex independently of interpolation.
            for(int a=0;a<3;++a) r[a]=sign(node,a);
            Point derivative;if(!evaluate(8,xyz,r,derivative)) return false;
            for(int k=0;k<24;++k) p.gradient[k]=derivative.gradient[k];
        }
    }
    p.left=node;return true;
}
MARS_MIXED_HD inline bool invert(int count,const double* xyz,const double* point,Point& p) {
    double r[3]={0,0,0};
    for(int iteration=0;iteration<20;++iteration) {
        double j[9],det;if(!evaluate(count,xyz,r,p,j,&det)) return false;
        double error[3],cofactor[9],largest=0;
        for(int a=0;a<3;++a) error[a]=point[a]-p.position[a];
        cross(j+3,j+6,cofactor);cross(j+6,j,cofactor+3);cross(j,j+3,cofactor+6);
        for(int a=0;a<3;++a) {double d=dot(error,cofactor+3*a)/det;r[a]+=d;largest=fmax(largest,fabs(d));}
        if(largest<1e-12) {
            if(!evaluate(count,xyz,r,p)) return false;
            for(int k=0;k<count;++k) if(p.shape[k]<-1e-10 || p.shape[k]>1+1e-10) return false;
            return true;
        }
    }
    return false;
}
}
#undef MARS_MIXED_HD
