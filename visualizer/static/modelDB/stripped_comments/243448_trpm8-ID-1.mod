NEURON {
    SUFFIX trpm8
	USEION ca READ ica, cai WRITE ica VALENCE 2
    RANGE minf, vhalf
    RANGE p_ca, accel
    RANGE em8, gbar, am8, C, z
	NONSPECIFIC_CURRENT im8
}


UNITS {
    R = (k-mole) (joule/degC)
    (mA) = (milliamp)
    (mV) = (millivolt)
    (mol) = (1)
    (molar) = (1/liter)
    (mM) = (millimolar)
} 

CONSTANT {
    F = 96500        (coulomb)        
}


PARAMETER {
	gbar = 1e-7        (mho/cm2)
    dE	 = 9e3           (joule)
    C	 = 67
    z	 = 0.65

    em8  = 0            (mV)
    mmin = 0		   (mV)
    mmax = 200		   (mV)
    Kca  = 0.0005       (mM)

    p_ca = 0.01    
    taum = 80000       (ms)
}

STATE {
    m  (mV)
}

INITIAL {
	rate(cai)
	m= minf
}

ASSIGNED {
    celsius (degC)
    v       (mV)
    ica     (mA/cm2)
    vhalf   (mV)
    am8
    minf    (mV)
	cai		(mM)
	im8     (mA/cm2)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	vhalf=(1000)*(C*R*celsius - dE)/(z*F)+m 
	am8=1/(1+exp(-z*F*(v-vhalf)/((1000)*R*(celsius+273.15)))) 
    im8 = (1-p_ca)*gbar*am8*(v-em8) 
	ica = p_ca*gbar*am8*(v-em8) 
	
}

DERIVATIVE states {
	rate(cai)
	m' = (minf-m)/taum
}

PROCEDURE rate(ca (mM)){
	minf = mmin+(mmax-mmin)*(ca)/(Kca+ca) 
}