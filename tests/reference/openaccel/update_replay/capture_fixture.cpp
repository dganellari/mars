#include "mars_update_reference.hpp"
#include <iostream>
int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);
    for(int iteration=0;iteration<2;++iteration) {
        mars_reference::UpdateSession session(MPI_COMM_WORLD);
        for(int phase=0;phase<8;++phase) {
            session.phase(phase);
            const int phases[]={2,2,5,3,3,3,4,1,1,1};
            // The mean must be available before trace application, as in the real routine.
            const int stages[]={0,1,2,3,4,5,6,7,9,8};
            for(int stage:stages) if(phases[stage]==phase) {
                const int samples=stage==3?6:(stage==4||stage==5||stage==7||stage==8?3:1);
                for(int sample=0;sample<samples;++sample) {
                    update_replay::Input in;in.stage=stage;double* x=in.values;
                    if(stage==0){x[0]=2;x[1]=1;x[2]=.3;}
                    if(stage==1)x[0]=1;
                    if(stage==2){x[0]=3;x[3]=2;x[4]=3;x[5]=4;x[6]=.1;x[7]=.2;x[8]=.3;}
                    if(stage==3||stage==5){x[0]=1;x[1]=3;x[19]=1.0/3;x[22]=1;x[23]=.75;}
                    if(stage==4){x[0]=1;x[1]=-3;x[4]=1.0/3;x[7]=-1;x[8]=.75;}
                    if(stage==6)for(int s=0;s<3;++s){x[s]=1;x[6+3*s]=3;x[15+s]=x[18+s]=2;x[21+3*s]=1.0/3;}
                    if(stage==7){for(int s=0;s<3;++s){x[s]=2;x[3+s]=1.0/3;}x[6]=1.0/3;}
                    if(stage==8){x[0]=x[1]=x[2]=x[4]=2;x[3]=.05;}
                    if(stage==9){x[0]=2;x[1]=1;}
                    update_replay::Output out;
                    update_reference::evaluate(in,out,stage==9?1:(stage==3?50:stage==4?20:stage>=5?30:10),sample);
                }
            }
        }
    }
    MPI_Finalize();std::cout<<"PASS: extracted capture expressions wrote two synthetic update sequences\n";
}
