UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
	SUFFIX hcurrent
	NONSPECIFIC_CURRENT i
	RANGE g, v50
	GLOBAL gfactor, eh
}
 
PARAMETER {
        v		(mV)
        celsius		(degC)
        g= 0.0001		(mho/cm2)
		
		v50=-82		(mV)
		gfactor = 1
}
 
STATE {
h
}
 
ASSIGNED {
        eh (mV)
	i		(mA/cm2) 
 	hinf
	htau    (ms)
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        i=g*h*(v-eh)*gfactor
}
 
DERIVATIVE states { 
       rates(v)
       h'= (hinf- h)/ htau
}

INITIAL { 
	rates(v)
	h= hinf
}


PROCEDURE rates(v (mV)) {
UNITSOFF





hinf = 1/(1+exp((v-v50)/10.5))
htau = (1/(exp(-14.59-0.086*v)+exp(-1.87+0.0701*v))) 
UNITSON
}