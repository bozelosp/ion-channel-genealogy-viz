UNITS {

	(mV)        = (millivolt) 
	(mA)        = (milliamp) 
	(pS) 		= (picosiemens)
	(um) 		= (micron)
		

} 

NEURON { 
	SUFFIX NaPyr
	USEION na READ ena WRITE ina
	RANGE  gbar, ina, minf,malpha, mbeta
	GLOBAL  v_table_min, v_table_max
}

PARAMETER {

	gbar =  	350     (pS/um2)
	v                   (mV)
	ena 		        (mV)  
	phih =  5	        (1)
	phim =	5           (1)
	

	v_table_min 		= -120			(mV)
	v_table_max 		= 100			(mV)
} 

ASSIGNED { 

	ina 		        (mA/cm2) 
	
	halpha              (1/ms)
	hbeta               (1/ms)
    malpha          	(1/ms)      
   	mbeta           	(1/ms)      
} 


STATE {
	h                 (1)
	m                 (1)
}

BREAKPOINT { 
	settables(v) 
	SOLVE states METHOD cnexp
	
    ina =(1e-4)*gbar * m * m* m * h * (v - ena)            
} 

INITIAL {
	settables(v) 
	h=halpha/(halpha+hbeta) 
	m=malpha/(malpha+mbeta)
} 

DERIVATIVE states { 
	h' =phih* ( halpha*(1-h) -hbeta*h) 
    m' =phim* ( malpha*(1-m) - mbeta*m)             
}


UNITSOFF


PROCEDURE settables(v (mV)) { 
	TABLE halpha, hbeta, malpha, mbeta FROM 	v_table_min  TO v_table_max WITH 961   
	
	halpha  = 0.128*(exp(-(v+50)/18))
	hbeta   = 4/(1+exp(-0.2*(v+27)))
	
	
	malpha = 0.32*vtrap(-(v+54),0.25)   
	mbeta =  0.28*vtrap((v+27),0.2)                      
	
	
}

UNITSON


FUNCTION vtrap(x, k) {
  if (fabs(x) < 1e-6) {
    vtrap = 1/(k * exp(k*x))
  } else {
    vtrap = x / (exp(k*x) - 1)
  }
}