
COMMENT

na12.mod

Sodium channel, Hodgkin-Huxley style kinetics.  

Kinetics were fit to data from Filipis et al. 2022

Author: Zach Mainen, Salk Institute, 1994, zach@salk.edu

ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX na16
	USEION na READ ena WRITE ina
	USEION ca READ eca WRITE ica
	USEION cal WRITE ical VALENCE 2
	RANGE m, h, gna, gbar,ina, ica, ical
	GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf, cana
	RANGE minf, hinf, mtau, htau
	GLOBAL Ra, Rb, Rd, Rg
	GLOBAL q10, temp, tadj, vmin, vmax, vshift, thi12, hw,hc
}
 
PARAMETER {
	cana=0
	gbar = 1000   	(pS/um2)	: 0.12 mho/cm2
	vshift = -2	(mV)		: voltage shift (affects all)
								
	tha  = -44	(mV)		: v 1/2 for act		(-42)
	qa   = 6.5	(mV)		: act slope		
	Ra   = 2.1	(/ms)		: open (v)		
	Rb   = 0.17	(/ms)		: close (v)		

	thi1  = -80	(mV)		: -35mV v 1/2 for inact 
	thi12  = -35	(mV)		: v 1/2 for inact 	
	thi2  = -75	(mV)		: v 1/2 for inact 	
	qi   = 5	(mV)	        : inact tau slope
	thinf  = -72	(mV)		: inact inf slope	
	qinf  = 6.2	(mV)		: inact inf slope
	Rg   = 0.0091	(/ms)		: inact (v)	
	Rd   = 0.024	(/ms)		: inact recov (v) 

	temp = 23	(degC)		: original temp 
	q10  = 2.3			: temperature sensitivity

	v 		(mV)
	dt		(ms)
	celsius		(degC)
	vmin = -120	(mV)
	vmax = 100	(mV)
	hc=0.01
	hw=0.02
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ina 		(mA/cm2)
			ica 		(mA/cm2)
			ical 		(mA/cm2)
		eca  (mV)
	gna		(pS/um2)
	ena		(mV)
	minf 
			
	hinf
	hinf2
	mtau (ms)	
	htau (ms)
	htau2 (ms)
	tadj
}
 

STATE { m h h2}

INITIAL { 
	trates(v+vshift)
	m = minf
	h = hinf
	h2=hinf2
}

BREAKPOINT {
        SOLVE states
        gna = tadj*gbar*m*m*m*h+tadj*gbar*m*m*m*hw*h2
	ina = (1e-4) * gna * (v - ena)
	ica = (1e-4) * gna * (v - eca)*cana/100
	ical=ica
} 

LOCAL mexp, hexp ,hexp2

PROCEDURE states() {   :Computes state variables m, h, and n 
        trates(v+vshift)      :             at the current v and dt.
        m = m + mexp*(minf-m)
        h = h + hexp*(hinf-h)
        h2 = h2 + hc
        VERBATIM
        //return 0;
        ENDVERBATIM
}

PROCEDURE trates(v) {  
                      
        LOCAL tinc
        TABLE minf, mexp, hinf, hexp, hinf2, hexp2
	DEPEND dt, celsius, temp, Ra, Rb, Rd, Rg, tha, thi1, thi2, qa, qi, qinf
	
	FROM vmin TO vmax WITH 199

	rates(v): not consistently executed from here if usetable == 1

        tadj = q10^((celsius - temp)/10)
        tinc = -dt * tadj

        mexp = 1 - exp(tinc/mtau)
        hexp = 1 - exp(tinc/htau)
        hexp2 = 1 - exp(tinc/htau2)
}


PROCEDURE rates(vm) {  
        LOCAL  a, b

	a = trap0(vm,tha,Ra,qa)
	b = trap0(-vm,-tha,Rb,qa)
	mtau = 1/(a+b)
	minf = a*mtau

		:"h" inactivation 

	a = trap0(vm,thi1,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	htau = 1/(a+b)

	hinf = 1/(1+exp((vm-thinf)/qinf))

	a = trap0(vm,thi12,Rd,qi)
	b = trap0(-vm,-thi2,Rg,qi)
	hinf2=hinf
	htau2 = hc
}


FUNCTION trap0(v,th,a,q) {
	if (fabs(v/th) > 1e-6) {
	        trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
	} else {
	        trap0 = a * q
 	}
}	




