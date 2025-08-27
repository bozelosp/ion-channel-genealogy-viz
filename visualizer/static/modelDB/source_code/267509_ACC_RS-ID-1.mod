COMMENT
//****************************//
// Created by Alon Polsky 	//
//    apmega@yahoo.com		//
//		2010			//
//****************************//
Modified 2015 by Robert Egger
to include facilitation variable
as modeled by Varela et al. 1997
ENDCOMMENT

TITLE NMDA synapse 


NEURON {
	POINT_PROCESS ACC_RS
	NONSPECIFIC_CURRENT inmda,iampa
	RANGE gampamax,gnmdamax,inmda,iampa
	RANGE gnmda,gampa
	RANGE e,tau1,tau2,tau3,tau4
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
	(mM)    = (milli/liter)
 	(mA) = (milliamp)
	(um) = (micron)
}

PARAMETER {
	gnmdamax=1	(nS)
	gampamax=1	(nS)
	e= 0.0	(mV)
	tau1= 60    (ms): NMDA inactivation
	tau2=4	(ms)	: NMDA activation
	tau3=	2    (ms)	: AMPA inactivation 
	tau4=0.2	(ms)	: AMPA activation
	
	n=0.05 	(/mM)	  
	gama=0.05	(/mV)   
	dt 		(ms)
	v		(mV)
	}

ASSIGNED { 
	inmda		(nA)  
	iampa		(nA)  
	gnmda		(nS)
	gampa		(nS)

}
STATE {
	A 		(nS)
	B 		(nS)
	C 		(nS)
	D 		(nS)
}


INITIAL {
    gnmda=0 
    gampa=0 
	A=0
	B=0
	C=0
	D=0
	
	}    

BREAKPOINT {  
    
	LOCAL count
	SOLVE state METHOD cnexp
	gnmda=(A-B)/(1+n*exp(-gama*v) )
	gampa=(C-D)
	inmda =(1e-3)*gnmda*(v-e)
	iampa= (1e-3)*gampa*(v-e)
	
}

NET_RECEIVE(weight_ampa, weight_nmda) {
 
	INITIAL {
	  gampamax = weight_ampa
	  gnmdamax = weight_nmda
	}
	gampamax = weight_ampa
	gnmdamax = weight_nmda
	
	A = A+ gnmdamax
	B = B+ gnmdamax
	C = C+ gampamax
	D = D+ gampamax

	    	
	VERBATIM
	printf("***********\n");
	//printf("A = %.2f\n", A);
	//printf("B = %.2f\n", B);
	//printf("C = %.2f\n", C);
	//printf("D = %.2f\n", D);
	ENDVERBATIM
}
DERIVATIVE state {
	A'=-A/tau1
	B'=-B/tau2
	C'=-C/tau3
	D'=-D/tau4
}