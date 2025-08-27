NEURON {
	SUFFIX dCaAP
	RANGE w, vth
	RANGE K, tauA,tauB, A, B, D
	RANGE t_dCaAP
	RANGE sigma_diff,refract_period
	RANGE vrest 
	RANGE dcaap_count
	NONSPECIFIC_CURRENT i
}

UNITS {
    (mV) = (millivolt)
    (nA) = (nanoamp)
    (uS) = (microsiemens)
}

PARAMETER {
    w = 0 
    vth    (mV)   
	refract_period
	tauA
	tauB
	D
	sigma_diff
	K
	vrest
	t_dCaAP 
	dcaap_count
}

ASSIGNED {
    v (mV)
   	i (nA)
}

STATE{
	A
	B
}
INITIAL {
	A=0
	B=0



	K=0 

	t_dCaAP = -refract_period 
	dcaap_count = 0
}

BREAKPOINT {	
	SOLVE dCaAP METHOD cnexp
	i =  -(A - B) * w * K
}

DERIVATIVE dCaAP{
	
	
	
	
	
	A' =  A * (1 - A) / tauA
	B' =  B * (1 - B) / tauB
}


BEFORE BREAKPOINT {
	
	
	
	if (t >= t_dCaAP + refract_period && v > vth && w > 0) {
		t_dCaAP = t
		A=0.001
		B=0	
		
		K = exp( -(v - vth) / (vth - vrest) / D ) 
		if(K > 1){
			K = 1
		}
		dcaap_count = dcaap_count + 1
		printf("Fired %.0f dCaAP at %f ms with K = %f\n",dcaap_count, t_dCaAP,K)
	}
	
	if(B == 0 && t >= t_dCaAP + sigma_diff && w > 0){
		B=0.001
	}		
}