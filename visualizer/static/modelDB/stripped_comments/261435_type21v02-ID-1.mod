UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (mS) = (millisiemens)
}
 
NEURON {
	SUFFIX type21
	NONSPECIFIC_CURRENT i
	RANGE  ninit            
	                        
	                        
	RANGE  type21           
	                        
	                        
	
		RANGE gl,el,v12,sl  
		RANGE n0,sn,t0,st,v0,sg
	
	
		RANGE gk,ek,gna,ena,a,b,cap, mhalf, mslope, S
	
	
	THREADSAFE 

}
 
PARAMETER {
        v               (mV)
        type21=2        (1)
        gna= 120.       (mS/cm2)
        ena=  50.       (mV)
        gk =  36.       (mS/cm2)
        ek = -77.       (mV)
        gl              (mS/cm2)
        el              (mV)
        n0              (1)
        sn              (1)
        t0              (ms)
        st              (ms)
        v0              (mV)
        sg              (mV)
        v12             (mV)
        sl              (mV)
        a  =  0.906483183915
        b  = -1.10692947808
        mhalf = -40
        mslope = 9.5
        cap = 1.0
        ninit = 0.34
        S = 1.3
	
	
}
 
STATE {
   n
}

ASSIGNED {
        i       (mA/cm2) 
        minf
        ninf
        ntau    (ms)
}
 
BREAKPOINT {
	SOLVE states METHOD cnexp
	
	i = (1e-3)*( gna*minf*minf*minf*(a+n*b)*(v-ena)+gk*pow(n/S,4.0)*(v-ek)+gl*(v-el) )/cap
}
 
DERIVATIVE states { 
	rates(v)
	n'= (ninf- n)/ ntau 
}


INITIAL {
	if ( fabs(type21 - 1.) < 1e-6 ){
		
		gl  =   0.3  (mS/cm2)
		el  = -54.3  (mV)
		n0  =   0.35
		sn  =   1. - n0
		v12 = -40.   (mV)
		sl  =   4.   (mV)
		t0  =    .46 (ms)
		st  =   3.5  (ms)
		v0  = -60.5  (mV)
		sg  =  35.9  (mV)
		
	} 
	if ( fabs(type21 - 2.) < 1e-6 ){
		
		gl  =   0.1  (mS/cm2)
		el  = -39.   (mV)
		n0  =   0.28
		sn  =   1. - n0
		v12 = -44.5  (mV)
		sl  =   9.   (mV)
		t0  =    .5  (ms)
		st  =   5.   (ms)
		v0  = -60.   (mV)
		sg  =  30.   (mV)
		
	} 
	rates(v)
	if (ninit < 0 || ninit > 1){
		n = ninf
	} else {
		n = ninit
	}
}

PROCEDURE rates(v (mV)) {
UNITSOFF 
	
	minf =      1./(1.+exp(-(v-mhalf)/mslope))
	ninf = n0 + sn/(1.+exp(-(v-v12)/sl ))
	
	ntau = t0/(exp((v-v0)/sg)+exp(-(v-v0)/sg)) 
UNITSON
}