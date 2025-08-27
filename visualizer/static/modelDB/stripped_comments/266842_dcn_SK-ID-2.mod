NEURON {
	SUFFIX dcnSK
	USEION ca READ cai VALENCE 2
	USEION k READ ek WRITE ik
	RANGE gbar, z, ik
	GLOBAL qdeltat
}
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}
 
PARAMETER { 
    qdeltat = 1
    gbar = 1e-5 (siemens/cm2)
} 

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2) 
	cai (mM)
	zinf 
    tauz (ms) 
} 
 
STATE {
	z 
} 

INITIAL { 
    rate(cai)
    z = zinf 
} 
 
BREAKPOINT { 
    SOLVE states METHOD cnexp 
	ik = gbar * z * (v - ek)
} 

DERIVATIVE states { 
	rate(cai) 
	z' = (zinf - z) / tauz
} 

PROCEDURE rate(cai(mM)) {
	TABLE zinf, tauz FROM 0 TO 0.01 WITH 300
    zinf = cai*cai*cai*cai / (cai*cai*cai*cai + 8.1e-15) 

    if (cai < 0.005) {
        tauz = 1 - (186.67 * cai)
    } else {
        tauz = 0.0667
    }
    tauz = tauz / qdeltat
}