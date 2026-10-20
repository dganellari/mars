// Included after the extracted device functions by test_amd_mfma_direct.py.
// The emulated MFMA uses the measured gfx942 lane layout; the reference below
// independently evaluates D^T G D with scalar tensor contractions.
std::vector<double> reference(const std::vector<double>& u, const std::vector<int>& dof,
                              const std::vector<double>& D, const std::vector<double>& G,
                              const std::vector<double>& weights, int elements, bool affine)
{
    std::vector<double> y(u.size(), 0);
    for(int e=0; e<elements; ++e) {
        double g[3][512] = {};
        for(int i=0;i<8;++i) for(int j=0;j<8;++j) for(int k=0;k<8;++k) {
            const int p=i*64+j*8+k;
            for(int s=0;s<8;++s) {
                g[0][p]+=D[i*8+s]*u[dof[e*512+s*64+j*8+k]];
                g[1][p]+=D[j*8+s]*u[dof[e*512+i*64+s*8+k]];
                g[2][p]+=D[k*8+s]*u[dof[e*512+i*64+j*8+s]];
            }
            const auto* m=G.data()+(affine?e*6:(e*512+p)*6);
            const double w=affine?weights[i]*weights[j]*weights[k]:1;
            const double a=g[0][p],b=g[1][p],c=g[2][p];
            g[0][p]=w*(m[0]*a+m[1]*b+m[2]*c);
            g[1][p]=w*(m[1]*a+m[3]*b+m[4]*c);
            g[2][p]=w*(m[2]*a+m[4]*b+m[5]*c);
        }
        for(int i=0;i<8;++i) for(int j=0;j<8;++j) for(int k=0;k<8;++k) {
            double v=0;
            for(int s=0;s<8;++s)
                v+=D[s*8+i]*g[0][s*64+j*8+k]+D[s*8+j]*g[1][i*64+s*8+k]+D[s*8+k]*g[2][i*64+j*8+s];
            y[dof[e*512+i*64+j*8+k]]+=v;
        }
    }
    return y;
}
template<int Mode>
bool check(int elements, bool constant)
{
    const int count=elements*256+256;
    std::vector<int> dof(elements*512);
    for(int e=0;e<elements;++e) for(int p=0;p<512;++p) dof[e*512+p]=e*256+p;
    std::vector<double> u(count),D(64),G(elements*(Mode?6:512*6)),weights(8),y(count+1,0);
    y.back()=12345;
    for(int p=0;p<count;++p) u[p]=constant?1:std::sin(.13*p)+.3*std::cos(.07*p);
    for(int i=0;i<8;++i) {
        weights[i]=.2+.1*i;
        double sum=0;
        for(int j=0;j<8;++j) if(i!=j) { D[i*8+j]=std::sin(1+i*8+j); sum+=D[i*8+j]; }
        D[i*8+i]=-sum;
    }
    for(size_t p=0;p<G.size()/6;++p) {
        const double scale=1+.1*std::sin(double(p));
        const double m[6]={2,.13,.07,3,.17,4};
        for(int j=0;j<6;++j) G[p*6+j]=scale*m[j];
    }
    const auto expected=reference(u,dof,D,G,weights,elements,Mode);
    for(int block=0;block<(elements+1)/2;++block) {
        blockIdx.x=block;
        std::vector<std::thread> threads;
        for(int lane=0;lane<64;++lane) threads.emplace_back([&,lane]{
            threadIdx.x=lane;
            mars::fem::ho_laplacian_apply_mfma_direct<Mode>(u.data(),y.data(),dof.data(),D.data(),G.data(),weights.data(),elements);
        });
        for(auto& t:threads) t.join();
    }
    double error=0,scale=1;
    for(int p=0;p<count;++p) {
        if(!std::isfinite(y[p])) return false;
        error=std::max(error,std::abs(y[p]-expected[p]));
        scale=std::max(scale,std::abs(expected[p]));
        if(constant && std::abs(y[p])>1e-10) return false;
    }
    const bool ok=error/scale<1e-12 && y.back()==12345;
    std::printf("%s lane emulator elements=%d metric=%d constant=%d relative=%.3e\n",ok?"PASS":"FAIL",elements,Mode,constant,error/scale);
    return ok;
}
int main() {
    bool ok=true;
    for(int n:{1,2,3}) for(bool constant:{false,true}) {
        ok=check<0>(n,constant)&&ok;
        ok=check<1>(n,constant)&&ok;
    }
    return ok?0:1;
}
