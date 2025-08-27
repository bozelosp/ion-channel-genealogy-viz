NEURON {
	SUFFIX iC
	USEION k READ ko, ki WRITE ik		
	USEION ca READ cai   
	RANGE  gkcbar,ik
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(S)  	= (siemens)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER {
        gkcbar	= 1.0	 (S/cm2)
	
}

ASSIGNED {
        v       (mV)
        cai	(mM)     
	celsius		 (degC)
	ek		(mV)
	ik	(mA/cm2)
	k1	(/ms)
	k2	(/ms)
	k3	(/ms)
	k4	(/ms)
	q10	(1)
	ko	(mM)
	ki	(mM)
}

STATE { cst ost ist }

BREAKPOINT { 
	SOLVE kin METHOD sparse
	ek=25*log(ko/ki)		
	ik = gkcbar * ost *( v - ek ) 
}

INITIAL {
	SOLVE kin STEADYSTATE sparse
}

KINETIC kin {
	rates(v, cai)
	~cst<->ost  (k3,k4)
	~ost<->ist  (k1,0.0)
	~ist<->cst  (k2,0.0)
	CONSERVE cst+ost+ist=1
}



PROCEDURE rates( v(mV), cai(mM)) {

	 k1=alp( 0.01, v,  -10.0,   1.0 ) 
	 k2=alp( 0.1, v, -120.0, -10.0 ) 

	 k3=alpha( 0.001, 1.0, v, -20.0, 7.0 ) *1.0e8* (cai*1.0(/mM) )^3
	 
	 k4=alp( 0.2, v, -44.0,  -5.0 ) 

}

FUNCTION alpha( tmin(ms), tmax(ms), v(mV), vhalf(mV), k(mV) )(/ms){
        alpha = 1.0 / ( tmin + 1.0 / ( 1.0 / (tmax-tmin) + exp((v-vhalf)/k)*1.0(/ms) ) )
}

FUNCTION alp( tmin(ms), v(mV), vhalf(mV), k(mV) )(/ms){
        alp = 1.0 / ( tmin + exp( -(v-vhalf) / k )*1.0(ms) )
}