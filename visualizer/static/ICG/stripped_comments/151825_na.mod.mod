INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na
	USEION na READ ena WRITE ina
	RANGE m, h, gna, gbar,vshiftm,vshifth,taum_scale,tauh_scale
	RANGE minf, hinf, mtau, htau
	GLOBAL a1,a2,a3,a4,v05m,km
	GLOBAL i1,i2,i3,i4,v05h,kh
	GLOBAL q10, temp, tadj, vmin, vmax
}

PARAMETER {
	gbar = 1.0   	(pS/um2)	
	vshiftm =-5	(mV)		
	vshifth =-10  (mV)		
	taum_scale= 1
	tauh_scale=1
								
	a1=0.058		(ms)		
	a2=0.114		(ms)
	a3=-36		(mV)
	a4=28			(mV)
	v05m=-38		(mV)
	km=10			(mV)

	i1=0.28		(ms)		
	i2=16.7 		(ms)
	i3=-60		(mV)
	i4=25			(mV)
	v05h=-66		(mV)
	kh=-6			(mV)

	
	temp = 21	(degC)		
	q10  = 2.3				

	v 		(mV)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
	gna		(pS/um2)
	ena		(mV)
	minf 		hinf
	mtau (ms)	htau (ms)
	tadj
}
 

STATE { m h }

INITIAL { 
	mrates(v+vshiftm)
	hrates(v+vshifth)
	m = minf
	h = hinf
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gbar*m*m*m*h
	ina = (1e-4) * gna * (v - ena)
} 


DERIVATIVE states {    

	mrates(v+vshiftm)
	hrates(v+vshifth)
        m' = (minf-m)/mtau
        h' = (hinf-h)/htau

}



PROCEDURE mrates(vm) {  

	

	tadj = q10^((celsius - temp)/10)
	mtau = (a1+a2*exp(-((vm-a3)/a4)^2))/tadj
	minf = 1/(1+exp(-(vm-v05m)/km))
}


PROCEDURE hrates(vm) {

	

        tadj = q10^((celsius - temp)/10)
	htau = (i1+i2*exp(-((vm-i3)/i4)^2))/tadj
	hinf = 1/(1+exp(-(vm-v05h)/kh))
}