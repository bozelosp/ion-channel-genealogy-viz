NEURON {
    SUFFIX na8stX
    USEION na READ ena WRITE ina
    GLOBAL vShiftX, vShift_inactX, maxrateX
    RANGE vShift_inact_local
    RANGE g, gbar, kinfact
    RANGE a1_0, a1_1, b1_0, b1_1, a2_0, a2_1
    RANGE b2_0, b2_1, a3_0, a3_1, b3_0, b3_1
    RANGE bh_0, bh_1, bh_2, ah_0, ah_1, ah_2
}

UNITS { (mV) = (millivolt) }



PARAMETER {
    gbar = 33     (millimho/cm2)

    a1_0 = 6.264774039489168e+01 (/ms)
    a1_1 = 1.160554780103536e-02 (/mV) 
    
    b1_0 = 1.936911472259165e-03 (/ms)
    b1_1 = 1.377185203515948e-01 (/mV)

    a2_0 = 3.478282276988217e+01 (/ms)
    a2_1 = 2.995594783341219e-02 (/mV) 
    
    b2_0 = 9.575149443481501e-02 (/ms)
    b2_1 = 9.281138012170398e-02 (/mV)

    a3_0 = 7.669829640279345e+01 (/ms)
    a3_1 = 5.374324331056838e-02 (/mV) 
    
    b3_0 = 1.248791525464647e+00 (/ms)
    b3_1 = 3.115037791363419e-02 (/mV)

    bh_0 = 3.573645069880386e+00 (/ms) 
    bh_1 = 1.933213300303968e-01 
    bh_2 = 7.496541077890667e-02 (/mV)

    ah_0 = 6.882666625638676e+00 (/ms)
    ah_1 = 4.654019001523467e+03  
    ah_2 = 2.958332680760088e-02 (/mV)

    vShiftX = 12            (mV)  
                                 
    vShift_inactX = 10      (mV)  
                                 
    vShift_inact_local = 0 (mV)  
    maxrateX = 8.00e+03     (/ms) 
                                 
    kinfact = 1               
}

ASSIGNED {
    v    (mV)
    ena  (mV)
    g    (millimho/cm2)
    ina  (milliamp/cm2)
    a1   (/ms)
    b1   (/ms)
    a2   (/ms)
    b2   (/ms)
    a3   (/ms)
    b3   (/ms)
    ah   (/ms)
    bh   (/ms)
}

STATE { c1 c2 c3 i1 i2 i3 i4 o }

BREAKPOINT {
    SOLVE kin METHOD sparse
    g = gbar*o
    ina = g*(v - ena)*(1e-3)
}

INITIAL { SOLVE kin STEADYSTATE sparse }

KINETIC kin {
    rates(v)
    ~ c1 <-> c2 (a1, b1)
    ~ c2 <-> c3 (a2, b2)
    ~ c3 <-> o (a3, b3)
    ~ i1 <-> i2 (a1, b1)
    ~ i2 <-> i3 (a2, b2)
    ~ i3 <-> i4 (a3, b3)
    ~ i1 <-> c1 (ah, bh)
    ~ i2 <-> c2 (ah, bh)
    ~ i3 <-> c3 (ah, bh)
    ~ i4 <-> o  (ah, bh)
    CONSERVE c1 + c2 + c3 + i1 + i2 + i3 + i4 + o = 1
}




PROCEDURE rates(v(millivolt)) {
    LOCAL vS
    vS = v-vShiftX

    a1 = kinfact*a1_0*exp( a1_1*vS)
    a1 = a1*maxrateX / (a1+maxrateX)
    b1 = kinfact*b1_0*exp(-b1_1*vS)
    b1 = b1*maxrateX / (b1+maxrateX)
    
    a2 = kinfact*a2_0*exp( a2_1*vS)
    a2 = a2*maxrateX / (a2+maxrateX)
    b2 = kinfact*b2_0*exp(-b2_1*vS)
    b2 = b2*maxrateX / (b2+maxrateX)
    
    a3 = kinfact*a3_0*exp( a3_1*vS)
    a3 = a3*maxrateX / (a3+maxrateX)
    b3 = kinfact*b3_0*exp(-b3_1*vS)
    b3 = b3*maxrateX / (b3+maxrateX)
    
    bh = kinfact*bh_0/(1+bh_1*exp(-bh_2*(vS-vShift_inactX-vShift_inact_local)))
    bh = bh*maxrateX / (bh+maxrateX)
    ah = kinfact*ah_0/(1+ah_1*exp( ah_2*(vS-vShift_inactX-vShift_inact_local)))
    ah = ah*maxrateX / (ah+maxrateX)
}