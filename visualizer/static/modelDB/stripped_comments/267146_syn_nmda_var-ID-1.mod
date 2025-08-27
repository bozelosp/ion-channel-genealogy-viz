NEURON {
    POINT_PROCESS syn_nmda_var
    RANGE tau_o                           
    RANGE tau_c                          
    RANGE erev                              
    RANGE syn_step
    RANGE i,var                                 
    NONSPECIFIC_CURRENT i
    THREADSAFE
    POINTER randObjPtr
}

UNITS {
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (S) = (siemens)
}

PARAMETER {
    tau_o = 5.0 (ms)
    tau_c = 80.0 (ms)
    erev = 0.0 (mV)
    c1 = 0.05
    c2 = -0.08
    syn_step = 1.25
}


ASSIGNED {
    v (mV)
    i (nA)
    var
    randObjPtr
}

STATE {
    o
    c
}

INITIAL {
    o = 0
    c = 0
}

PROCEDURE seed(x) {
    set_seed(x)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    i = mgBlock_std(v) * (c - o) * (v-erev)
}

NET_RECEIVE(weight (uS)) {
    var = ceil(5*randGen())
    
    
    o = o + syn_step*weight*2*(var-1)/10.0
    c = c + syn_step*weight*2*(var-1)/10.0
}

VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION randGen() {
VERBATIM
	if (_p_randObjPtr) {
		
		_lrandGen = nrn_random_pick(_p_randObjPtr);
	}else{
		hoc_execerror("multithread random in NetStim"," only via hoc Random");
	}
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
	void** pv = (void**)(&_p_randObjPtr);
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
ENDVERBATIM
}

DERIVATIVE states {
    o' = -o/tau_o
    c' = -c/tau_c
}

UNITSOFF
FUNCTION mgBlock_std(v (mV)) {
    mgBlock_std = 1 / (1 + c1*exp(c2*v))
}
UNITSON