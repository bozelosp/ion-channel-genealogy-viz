NEURON {
    SUFFIX na8st
    USEION na READ ena WRITE ina
    GLOBAL vShift, vShift_inact, slow
    RANGE vShift_inact_local
    RANGE g, gbar
    RANGE a1_0, a1_1, b1_0, b1_1, a2_0, a2_1
    RANGE b2_0, b2_1, a3_0, a3_1, b3_0, b3_1
    RANGE bh_0, bh_1, bh_2, ah_0, ah_1, ah_2
}

UNITS { (mV) = (millivolt) 
(S) = (siemens)
}



PARAMETER {
    gbar = 0     (S/cm2)
	slow = 0
    a1_0 =0 (/ms) 
    a1_1 = 0 (/mV) 
    
    b1_0 = 0 (/ms) 
    b1_1 = 0 (/mV) 

    a2_0 = 0 (/ms) 
    a2_1 = 0 (/mV) 
    
    b2_0 = 0 (/ms) 
    b2_1 = 0 (/mV) 

    a3_0 = 0 (/ms) 
    a3_1 = 0 (/mV) 
    
    b3_0 = 0 (/ms) 
    b3_1 = 0 (/mV) 

    bh_0 = 0 (/ms) 
    bh_1 = 0
    bh_2 = 0 (/mV) 

    ah_0 = 0 (/ms) 
    ah_1 = 0
    ah_2 =0 (/mV) 

    vShift = 12            (mV)  
                                 
    vShift_inact = 10      (mV)  
                                 
    vShift_inact_local = 0 (mV)  
    maxrate = 8.00e+03     (/ms) 
                                 
}

ASSIGNED {
    v    (mV)
    ena  (mV)
    g    (S/cm2)
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

STATE { c1 c2 c3 i1 i2 i3 i4 i5 i6 o }

BREAKPOINT {
    SOLVE kin METHOD sparse
    g = gbar*o
    ina = (g)*(v - ena)
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
	~ i5 <-> c3 (ah/10, slow*bh/10)
    ~ i6 <-> o  (ah/10, slow*bh/10)
    CONSERVE c1 + c2 + c3 + i1 + i2 + i3 + i4 + i5 + i6 + o = 1 
}


PROCEDURE rates(v(millivolt)) {
    LOCAL vS
    vS = v-vShift
	
    a1 = a1_0*exp( a1_1*vS)
    b1 = b1_0*exp(-b1_1*vS)

    
    a2 = a2_0*exp( a2_1*vS)
    b2 = b2_0*exp(-b2_1*vS)

    
    a3 = a3_0*exp( a3_1*vS)
    b3 = b3_0*exp(-b3_1*vS)

    
    bh = bh_0/(1+bh_1*exp(-bh_2*(vS-vShift_inact-vShift_inact_local)))

    ah = ah_0/(1+ah_1*exp( ah_2*(vS-vShift_inact-vShift_inact_local)))
	

		a1 = a1*maxrate / (a1+maxrate)
		b1 = b1*maxrate / (b1+maxrate)
		a2 = a2*maxrate / (a2+maxrate)
		b2 = b2*maxrate / (b2+maxrate)
		a3 = a3*maxrate / (a3+maxrate)
		b3 = b3*maxrate / (b3+maxrate)
		bh = bh*maxrate / (bh+maxrate)
		ah = ah*maxrate / (ah+maxrate)
}