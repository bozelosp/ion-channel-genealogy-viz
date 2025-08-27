UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

}

PARAMETER {
	v 		(mV)
	eh  		(mV)        
	celsius 	(degC)
	gmax=.0001 	(mho/cm2)
	vhalfl=-90   	(mV)  
	vhalft=-75   	(mV)  
	a0t= .007 (/ms) 
	zetal= 2 (1) 
	zetat= 1.1 (1) 
	gmt=.4 (1)
	q10=4.5
	qtl=1
	myslope=0.07 
}


NEURON {
	SUFFIX ch_HCNp
	NONSPECIFIC_CURRENT i
	RANGE gmax, vhalfl, myi
	GLOBAL linf,taul, eh
}

STATE {
	l
}

ASSIGNED {
	i (mA/cm2)
	myi (mA/cm2)
	linf      
	taul
	g
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gmax*l
	i = g*(v-eh)
	myi = i
}


FUNCTION alpl(v(mV)) {
	alpl = exp(myslope*zetal*(v-vhalfl)) 
}

FUNCTION alpt(v(mV)) {
	alpt = exp(myslope*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
	bett = exp(myslope*zetat*gmt*(v-vhalft)) 
}

DERIVATIVE states {     
	rate(v)
	l' =  (linf - l)/taul
}

PROCEDURE rate(v (mV)) { 
	LOCAL a,qt
	qt=q10^((celsius-33)/10)
	a = alpt(v)
	linf = 1/(1+ alpl(v))
	taul = bett(v)/(qtl*qt*a0t*(1+a))
}