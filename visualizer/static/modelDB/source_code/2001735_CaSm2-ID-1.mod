TITLE Motoneuron Soma channels
: Calcium channels + Calcium Dynamics - Soma
: Updated by Mohamed Hisham , to correct the calcium dynamics
: SK channel time constant for the ADP


NEURON {
	SUFFIX CaSm2

	NONSPECIFIC_CURRENT iSK3
	NONSPECIFIC_CURRENT iSK2
	NONSPECIFIC_CURRENT ican
	NONSPECIFIC_CURRENT ical

	RANGE gSK3bar, gSK2bar, gcanbar, gcalbar, eca , ek
	RANGE gSK3, gSK2 , gcan, gcal
	RANGE mn_inf, hn_inf, ml_inf , S3_inf , S2_inf
	:RANGE tau_mn, tau_hn, tau_ml
	RANGE tmn, thn , S3_tau , S2_tau
	RANGE nexp3, kd3 , nexp2, kd2
	RANGE shiftT
	RANGE f, kca, alpha
	RANGE thetamn, thetahn
}


UNITS {
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
	(molar) = (1/liter)
	(mM) 	= (millimolar)

	FARADAY	= (faraday) (coulomb)
	R		= (k-mole) (joule/degC)
}

PARAMETER {

	: Calcium N-type channels
	gcanbar = 0.072837  (mho/cm2)	
	tmn		= 15	    (ms)
	thetamn	= 22	    (mV)
	thn		= 50	    (ms)
	thetahn	= 40	    (mV)

	: Calcium L-type Channels
	gcalbar = 0.0002    (mho/cm2)	
	mlexp	= 1
	tml		= 400	    (ms)
	thetaml = 45.8	    (mV)
	kml		= -3.7	    (mV)

	: Calcium-activated Potassium Channels SK3
	gSK3bar = 0.37418   (mho/cm2)	                    
	nexp3	= 1
	kd3		= 0.0002
	S3_tau 	= 40		(ms)

	: Calcium-activated Potassium Channels SK2
	gSK2bar = 0.37418   (mho/cm2)	                    
	nexp2	= 1
	kd2		= 0.0002
	S2_tau 	= 40		(ms)

	: Calcium Dynamics
	cao		= 2	    	(mM)
	caio	= .0000001  (mM)
	f		= 0.025
	alpha	= 0.08
	kca		= 0.7

	: General
	vtraub		= -10	    	(mV)
	celsius 	= 36        	(degC)
	shiftT 		= 0				(degC)
	ek 			= -90           (mV)
	eca 		= 134 			(mV)
}

STATE {
	mn hn ml S3 cai S2
}

ASSIGNED {
	dt      (ms)
	v       (mV)
	:eca		(mV)
	
	ican	(mA/cm2)
	ical	(mA/cm2)
	iSK3	(mA/cm2)
	iSK2	(mA/cm2)

	gSK3	(mho/cm2)
	gSK2	(mho/cm2)
	gcan	(mho/cm2)
	gcal	(mho/cm2)
	
	mn_inf
	hn_inf
	ml_inf
: calcium activated Potassium channel sKl
	S3_inf
	S2_inf
	
	tau_mn	(ms)
	tau_hn	(ms)
	tau_ml	(ms)
	
	tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp

	:eca = ((1000 * R * (celsius + 273.15)) / (2 * FARADAY)) * log(cao/cai)

	gcan = gcanbar * mn*mn*hn
	ican = gcan * (v - eca)

	gcal = gcalbar * (ml^mlexp)
	ical = gcal * (v - eca)

	gSK3 = gSK3bar * S3
	iSK3 = gSK3 * (v - ek)

	gSK2 = gSK2bar * S2
	iSK2 = gSK2 * (v - ek)
}

DERIVATIVE states {   : exact Hodgkin-Huxley equations

    evaluate_fct(v)

	mn' = (mn_inf - mn) / tau_mn
	hn' = (hn_inf - hn) / tau_hn
	ml' = (ml_inf - ml) / tau_ml
	S3'	= (S3_inf  - S3 ) / S3_tau
	S2'	= (S2_inf  - S2 ) / S2_tau

	cai' = f*(-(alpha*(ican+ical))-(kca*cai))
}

UNITSOFF
INITIAL {

	:  Q10 was assumed to be 3
	tadj = 3.0 ^ ((celsius-36-shiftT)/ 10 )

	cai = caio

	evaluate_fct(v)

	mn 	= mn_inf
	hn 	= hn_inf
	ml 	= ml_inf
	S3  = S3_inf
	S2  = S2_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL v2

	v2 = v - vtraub : convert to traub convention
	
	tau_mn = tmn * tadj
	mn_inf = 1 / (1+exp((v2+thetamn)/-5))
	
	tau_hn = thn * tadj
	hn_inf = 1 / (1+exp((v2+thetahn)/5))
	
	tau_ml = tml * tadj
	ml_inf = 1 / (1+exp((v2+thetaml)/kml))

:SK dynamics
	:S3_inf = ( (cai^nexp3) / ((cai^nexp3)+kd3) ) : just correct at nexp3 = 1
	S3_inf = 1/( 1+ (kd3/cai)^nexp3 )
	: SK2.
	S2_inf = 1/( 1+ (kd2/cai)^nexp2 )

}

