#include "mars_update_reference.hpp"
#include <iostream>
#include <random>
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    std::mt19937 generator(41); std::uniform_real_distribution<double> random(-3,3);
    std::size_t checks=0;double worst=0;
    for(int trial=0;trial<300;++trial)for(int stage=0;stage<10;++stage) {
        update_replay::Input in;in.stage=stage;
        for(double& value:in.values)value=random(generator);
        double* x=in.values;
        if(stage==0)x[2]=.3;
        if(stage==2)x[9]=trial%2;
        if(stage==3 || stage==5){x[0]=2+std::abs(x[0]);x[23]=.75;x[24]=stage==5 ? trial%2:0;}
        if(stage==4){x[0]=1+std::abs(x[0]);x[8]=.75;}
        if(stage==6){x[3]=x[4]=x[5]=trial%2;x[30]=(trial/2)%2;}
        if(stage==7)x[9]=trial%2;
        if(stage==8){x[3]=.05;x[5]=trial%2;}
        if(stage==9)x[1]=1+std::abs(x[1]);
        update_replay::Output expected,actual;
        update_reference::evaluate(in,expected);update_replay::evaluate(in,actual);
        for(int j=0;j<6;++j) {
            ++checks;const double error=std::abs(actual.values[j]-expected.values[j])/std::max(1.0,std::abs(expected.values[j]));
            worst=std::max(worst,error);
            if(!std::isfinite(actual.values[j]) || error>1e-12){std::cerr<<"native mismatch stage="<<stage<<" entry="<<j<<'\n';MPI_Abort(MPI_COMM_WORLD,1);}
        }
    }
    std::cout<<"PASS: "<<checks<<" comparisons against pinned update expressions; worst_scaled="<<worst<<'\n';
    MPI_Finalize();
}
