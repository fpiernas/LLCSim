# LLCSim
C++ Simulation Library for LLC resonant converters by Francisco Piernas Díaz

This is a simple header library to perform time domain simulations for Half-Bridge and Full-Bridge resonant LLC converters.

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

    You can change frequency on simulation time by calling circuit1.set_freq(unsigned int freq).
    It is recommended to change frequency only when previous period ends
    
    
   To change the load call change_R(double R);
