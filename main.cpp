/*
    LLC Simulation program + Q-Learning control
    Fran Piernas Diaz
    Compile with -O2

    Q-Learner tested paramsters: Alpha=0.7, Epsilon=0.0, Gamma=0.15
    Epsilon must be changed to higher value after enough iterations
*/
#include <iostream>
#include <cstdlib>
#include <fstream>
#include "llc.h"
using namespace std;
int main()
{
    ofstream plotter;
    plotter.open("sim.dat");
    llc circuit1(13e-9,150e-6,448.512e-6,1.752e-6,600e-6,0.86,1e-9,385,LLC_HALF_BRIDGE,80e3);
    for(unsigned long int i=0;circuit1.state.t<=4e-3;i++)
    {
        circuit1.step();
        if(i%10000000==0) cout<<"t= "<<circuit1.state.t<<" Vo= "<<circuit1.Vout<<endl;
        if(i%50==0)plotter<<circuit1.state.t<<" "<<circuit1.Vout<<endl;
        if(i==2000000) circuit1.set_freq(80.5e3);
    }
    return 0;
}
