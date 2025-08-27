NEURON {
    POINT_PROCESS syn_ampa_var
    RANGE tau_o                           : parameter
    RANGE tau_c                          : parameter
    RANGE erev                              : parameter
    RANGE syn_step
    RANGE i,var                                 : exposure
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
    tau_o = 0.2 (ms)
    tau_c = 3.0 (ms)
    erev = 0.0 (mV)
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
    i = (c - o) * (v-erev)
}

NET_RECEIVE(weight (uS)) {
    var = ceil(5*randGen()) 
    :var = ceil(5*scop_random()) : not thread safe!!
    :printf("%g \n", 2*(var-1)/10.0)
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
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.uniform(0,1)
		*/
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
