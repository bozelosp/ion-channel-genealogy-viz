NEURON {
	SUFFIX CaDen

	NONSPECIFIC_CURRENT ikca
	NONSPECIFIC_CURRENT ical

	RANGE gkcabar, gcabar, eca , gkca 
	RANGE tRandG , AmpRandG, O_inf, Warm_Gear, Warm_thresh ,tailon ,W_tau_d
}


UNITS {
	(mA)		= (milliamp)
	(mV)		= (millivolt)
	(molar)	= (1/liter)
	(mM)		= (millimolar)

	FARADAY	= (faraday) (coulomb)
	R			= (k-mole) (joule/degC)
}

PARAMETER {

	
	gcabar		= 0.0003		(mho/cm2)
	eca			= 60			(mV) 
	theta_m		= -30     	(mV)
	kappa_m		= -6	   	(mV)
	O_tau  		= 20      	(ms)
	O_tau2 		= 50			(ms)   
	W_tau_d    	= 1200		(ms)
	tailon		= 1			
	AmpRandG	= 0.5		
	Warm_Gear	= 1			
	Warm_thresh	= 0.27


	
	gkcabar		= 0.37418		(mho/cm2)
	nexp		= 10				(1)       
	kd			= 0.0005			(mM)
	S_tau		= 40				(ms)

	
	caio		= 0.0001			(mM)   
	f			= 0.01
	alpha		= 1				(cm2 mM/mC)
	kca	  		= 8				(1/ms)


	
	celsius 	= 36				(degC)
	ek      	= -80				(mV)
}


STATE {
	O W S cai (mM)
}

ASSIGNED {
	dt      (ms)
	v       (mV)

	
	O_inf

	
	S_inf 		(ms)

	
	W_inf
	W_tau  		(ms)
	W_tau2 		(ms)

	
	tRandG				
	

	
	ical	(mA/cm2)
	ikca	(mA/cm2)

	gkca				
	df 					
}

BEFORE BREAKPOINT {
	
	
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ical = gcabar * ((O)*(AmpRandG*0.3 + 1)) * (v - eca)   

	gkca = S
	ikca = gkcabar *  gkca * (v - ek)					
}

DERIVATIVE states {   

		rates(v)					
		S_inf 	= ( 1 / ((kd/(cai-0.00001))^nexp + 1) )
		df 		= d_Flag(O_inf-O)	
		tauCorrection(tRandG)		

		
		O'	=    (O_inf - O)*tauFunc(W,df) / (O_tau + O_tau2 * df) 
		W'	=    ( (W_inf*(1-df)) - W)/ (W_tau + W_tau2 )   
		S'	=	 (S_inf-S)/S_tau							
		cai' = f*(-(alpha*(ical))-(kca*cai))  
}

UNITSOFF

INITIAL {
		rates(v)
		O   = O_inf
		W   = W_inf
		cai = caio
		S   = 0
}

PROCEDURE rates(v(mV)) {
		TABLE O_inf , W_inf , W_tau
				FROM -200 TO 100 WITH 300

		O_inf = 1/(1+ exp((v - theta_m)/kappa_m) )
		W_inf = 1/(1+exp(-(v+57)/0.8))
		W_tau = 50 + (1150/Warm_Gear)/(1+exp((v+32)/7))
}



FUNCTION tauFunc(W,df){ LOCAL W_thresh, observable_W
		W_thresh			= Warm_thresh 
		observable_W	= tailon * W * df		
		
		
		tauFunc			= 0.001 + 0.999/(1+ exp((observable_W - W_thresh)/0.006)) 
}

PROCEDURE RandGenerator(O,O_inf){
	if( O<0.03 && O_inf > 0.025){
		tRandG   = scop_random() 
		AmpRandG = scop_random() 
		
	}
}


PROCEDURE tauCorrection(tRandG){ LOCAL tauCompansate
	 	tauCompansate		= (W_tau_d * (1 -  (0.76 *tRandG)) ) - W_tau
	 	W_tau2				= tauCompansate * df
}

FUNCTION d_Flag(Xf){
		
		d_Flag 	= 1/(1+expmod(5000*(Xf+0.01)))  
}

FUNCTION expmod(x){
	if(x<(-100)){
		expmod = 0
	}else{
		expmod = exp(x)
	}
}

UNITSON