/**
    LLC Simulation Library V3

    @author Piernas Diaz, Fran
    University of Granada, Spain

    MIT License

    Copyright (c) 2020 Francisco Piernas Diaz

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

    Eigen Library needed, compile with -O2

    Create a circuit with:

        llc circuit1(Cr,Lr,Lm,Ls,C,R,Ts,V,type,freq);

        where:
            Cr: resonant capacitor value, example 13e-9
            Lr: resonant inductor value, example 150e-6
            Lm: primary inductance value, example 448e-6
            Ls: secondary inductance value, example 1.7e-6
            C: output capacitor value, example 600e-6
            R: load in ohms, example 0.86
            Ts: time step, 1e-9 should be enough
            V: input voltage, example 385
            type: LLC_FULL_BRIDGE or LLC_HALF_BRIDGE
            freq: working frequency

    Once the object is created, call successively circuit1.step() to advance a time step.

    The circuit created is an equivalent version of the real circuit where C and R are
    reflected to primary.

    Check circuit status with circuit1.state.stateVector(i), with i being:
        0: Cr voltage
        1: Resonant current
        2: Primary current
        3: Output voltage

    Check real output voltage with circuit1.Vout

    Change frequency calling circuit1.set_freq(unsigned int freq).
    It is recommended to change frequency only when previous period ends


*/


#ifndef LLC_H_INCLUDED
#define LLC_H_INCLUDED

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <Eigen/Core>

using namespace std;
using namespace Eigen;

#define LLC_FULL_BRIDGE 1
#define LLC_HALF_BRIDGE 0
#define BRIDGE_BLOCKING 0
#define BRIDGE_FORWARD 1
#define BRIDGE_REVERSE -1
#define SLOPE 20e-6
class oscillator
{
    public:
        double t0;
        unsigned int freq;
        oscillator()
        {
            t0=0;
        }
        bool get_val(double& time)
        {
            double int_part;
            if(modf(1.0*freq*(time-t0),&int_part)<=0.5) return true;
            else return false;
        }
        void change_freq(unsigned int& freq, double now)
        {
            this->freq=freq;
            t0=now;
            return;
        }
};

double modv(double& V,double& time)
{
    if(time/SLOPE>1.0) return V;
    else return V/SLOPE*time;
}

struct vars
{
    Vector4d stateVector;
    double t=0; //simulation time
    bool pullupMosActivated=true;
    signed char bridge_state=BRIDGE_BLOCKING;
};

struct parameters
{
    double Cr=0;
    double Lr=0;
    double Lp=0;
    double Ls=0;
    double C=0;
    double Ceq=0;
    double Req=0;
    double R=0;
    double time_step=0;
    double V=0;
    bool converter_type=LLC_FULL_BRIDGE;
    unsigned int switching_freq=0;
    double M;
    double n; //turns ratio primary/secondary
};

class llc
{
    private:
        oscillator osc;
        Matrix4d A_forward, A_blocking, A_reverse;
        Vector4d B_forward, B_blocking, B_reverse;

        void update_State_Matrix(void)
        {
            A_forward=Matrix4d::Zero();
            A_blocking=Matrix4d::Zero();
            A_reverse=Matrix4d::Zero();

            B_forward=Vector4d::Zero();
            B_blocking=Vector4d::Zero();
            B_reverse=Vector4d::Zero();

            //Configure A_Forward:
            A_forward(0,1)=1.0/params.Cr;
            A_forward(1,0)=-1.0/params.Lr;
            A_forward(1,3)=-1.0/params.Lr;
            A_forward(2,3)=1.0/params.Lp;
            A_forward(3,1)=1.0/params.Ceq;
            A_forward(3,2)=-1.0/params.Ceq;
            A_forward(3,3)=-1.0/(params.Ceq*params.Req);

            //Configure A_blocking:
            A_blocking(0,1)=1.0/params.Cr;
            A_blocking(1,0)=-1.0/(params.Lr+params.Lp);
            A_blocking(2,0)=-1.0/(params.Lr+params.Lp);
            A_blocking(3,3)=-1.0/(params.Ceq*params.Req);

            //Configure A_reverse:
            A_reverse(0,1)=1.0/params.Cr;
            A_reverse(1,0)=-1.0/params.Lr;
            A_reverse(1,3)=1.0/params.Lr;
            A_reverse(2,3)=-1.0/params.Lp;
            A_reverse(3,1)=-1.0/params.Ceq;
            A_reverse(3,2)=1.0/params.Ceq;
            A_reverse(3,3)=-1.0/(params.Ceq*params.Req);

            //configure B_forward:
            B_forward(1)=1.0/params.Lr;

            //configure B_blocking:
            B_blocking(1)=1.0/(params.Lr+params.Lp);
            B_blocking(2)=1.0/(params.Lr+params.Lp);

            //configure B_reverse:
            B_reverse(1)=1.0/params.Lr;

            //Calculate equivalent dicrete matrix
            A_reverse=(Matrix4d::Identity()+params.time_step*A_reverse);
            A_blocking=(Matrix4d::Identity()+params.time_step*A_blocking);
            A_forward=(Matrix4d::Identity()+params.time_step*A_forward);

            B_forward=B_forward*params.time_step;
            B_blocking=B_blocking*params.time_step;
            B_reverse=B_reverse*params.time_step;

            return;
        }

    public:
        vars state; //state of equivalent circuit
        parameters params;
        //Variables for Vout of real circuit (not the equivalent simulate one)
        double Vout=0, Vout_prev=0;
        double Vout_der=0;

        void change_R(double R)
        {
            params.R=R;
            params.Req=params.R*params.n*params.n;
            update_State_Matrix();
            return;
        }


        llc(double Cr, double Lr, double Lp, double Ls, double C, double R, double time_step, double V, bool converter_type, unsigned int switching_freq)
        {
            params.Cr=Cr;
            params.Lr=Lr;
            params.Lp=Lp;
            params.Ls=Ls;
            params.C=C;
            params.R=R;
            params.time_step=time_step;
            params.V=V;
            params.converter_type=converter_type;
            params.switching_freq=switching_freq;
            params.M=sqrt(Lp*Ls);
            params.n=sqrt(Lp/Ls);
            params.Ceq=params.C/(params.n*params.n);
            params.Req=params.R*params.n*params.n;
            osc.change_freq(switching_freq,state.t);
            state.stateVector=Vector4d::Zero();
            update_State_Matrix();

            return;
        }

        void set_freq(unsigned int freq)
        {
            params.switching_freq=freq;
            osc.change_freq(freq,state.t);
            return;
        }

        void step(void)
        {
            state.pullupMosActivated=osc.get_val(state.t);
            double vg;
            if(params.converter_type==LLC_FULL_BRIDGE)
            {
                //Full bridge
                if(state.pullupMosActivated) vg=modv(params.V,state.t);
                else vg=-modv(params.V,state.t);
            }
            else
            {
                //Half bridge
                if(state.pullupMosActivated) vg=modv(params.V,state.t);
                else vg=0;
            }
            //StateVector= (Vcr, Ir, Ip, Vo)
            double ip_prev=state.stateVector(2);
            if(state.bridge_state==BRIDGE_FORWARD) state.stateVector=A_forward*state.stateVector+B_forward*vg;
            else if(state.bridge_state==BRIDGE_BLOCKING) state.stateVector=A_blocking*state.stateVector+B_blocking*vg;
            else if(state.bridge_state==BRIDGE_REVERSE) state.stateVector=A_reverse*state.stateVector+B_reverse*vg;

            double lp_v=(state.stateVector(2)-ip_prev)/params.time_step*params.Lp;
            if(state.bridge_state==BRIDGE_FORWARD&&(state.stateVector(1)-state.stateVector(2))<=0)
            {
                state.stateVector(1)=state.stateVector(2);
                state.bridge_state=BRIDGE_BLOCKING;
            }
            else if(state.bridge_state==BRIDGE_BLOCKING&&lp_v>state.stateVector(3)) state.bridge_state=BRIDGE_FORWARD;
            else if(state.bridge_state==BRIDGE_BLOCKING&&lp_v<-state.stateVector(3))state.bridge_state=BRIDGE_REVERSE;
            else if(state.bridge_state==BRIDGE_REVERSE&&(state.stateVector(1)-state.stateVector(2))>=0)
            {
                state.stateVector(1)=state.stateVector(2);
                state.bridge_state=BRIDGE_BLOCKING;
            }

            Vout=state.stateVector(3)/params.n;
            Vout_der=(Vout-Vout_prev)/params.time_step;
            Vout_prev=Vout;
            state.t+=params.time_step;
            return;
        }

        void next(void)
        {
            if(state.pullupMosActivated==true)
            {
                while(state.pullupMosActivated==true) step();
                while(state.pullupMosActivated==false) step();
            }
            else while(state.pullupMosActivated==false) step();
            return;
        }
};


#endif // LLC_H_INCLUDED
