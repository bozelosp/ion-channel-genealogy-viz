UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S)  = (siemens)
    (pS) = (picosiemens)
    (um) = (micrometer)
} 


NEURON {
    SUFFIX HH
    USEION na READ ena WRITE ina
    USEION k READ ek WRITE ik


    RANGE gnabar, gkbar
    RANGE gna, gk
	RANGE nm,  gkmbar
	RANGE m_exp, h_exp, n_exp, km_exp
	RANGE km_inf, tau_km
    RANGE m_inf, h_inf, n_inf	

    GLOBAL seed    

	GLOBAL vshift

	GLOBAL taukm,NF
    THREADSAFE
} 


PARAMETER {
    gnabar  = 0.12   (S/cm2)
    gkbar   = 0.036  (S/cm2)
	gkmbar = .002   	(S/cm2)
    
	
	
    gamma_na = 10  (pS)		
    gamma_k  = 10  (pS)		
    gamma_km  = 10  (pS)		
    seed = 1              
	vshift=0				
	taukm=1					
	NF=1					
	hslow=100
	hfast=0.3

} 


STATE {
    m h n nm 

} 


ASSIGNED {
    ina        (mA/cm2)
	ikdr	(mA/cm2)
	ikm	 (mA/cm2)
	ik		(mA/cm2)
	
    gna        (S/cm2)
    gk         (S/cm2)
	gkm         (S/cm2)
    ena        (mV)
    ek         (mV)
   
    dt         (ms)
    area       (um2)
    celsius    (degC)
    v          (mV)
        
    m_exp h_exp n_exp km_exp 
    m_inf h_inf n_inf km_inf 
    tau_m (ms) tau_h (ms) tau_n (ms) tau_km (ms) 



} 

INITIAL {
    
    rates(v)
    m = m_inf
    h = h_inf
    n = n_inf
	nm=km_inf
	

    set_seed(seed)
} 


BREAKPOINT {
    SOLVE states
    gna = gnabar * ((m)*(m)*(m)*(h))

    ina = gna * (v - ena)
    gk = gkbar * ((n)*(n)*(n)*(n))

    ikdr  = gk * (v - ek)
	gkm=gkmbar*((nm))

	ikm = gkm*(v-ek)
	ik=ikm+ikdr
	

	
} 


PROCEDURE states() {
    rates(v+vshift)
	m = m + m_exp * (m_inf - m)
	h = h + h_exp * (h_inf - h)
    n = n + n_exp * (n_inf - n)
	nm = nm + km_exp*(km_inf-nm)


    VERBATIM
    return 0;
    ENDVERBATIM
} 


PROCEDURE rates(vm (mV)) { 
    LOCAL a,b,i
    
    UNITSOFF
    
    
 
	a =-.6 * vtrap((vm+30),-10)	
	b = 20 * (exp((-1*(vm+55))/18))
	tau_m = 1 / (a + b)
	m_inf =addnoise( a * tau_m)

 

	a = 0.4 * (exp((-1*(vm+50))/20))
	b = 6 / ( 1 + exp(-0.1 *(vm+20)))
	tau_h=hslow/((1+exp((vm+30)/4))+(exp(-(vm+50)/2)))+hfast
	h_inf=addnoise( 1/(1 + exp((vm + 44)/4)))
	
   

	a = -0.02 * vtrap((vm+40),-10)
	b = 0.4 * (exp((-1*(vm + 50))/80))
	tau_n = 1 / (a + b)
	n_inf =addnoise( a * tau_n)


    a = -.001/taukm * vtrap((vm+30),-9)
    b =.001/taukm * vtrap((vm+30),9) 
    tau_km = 1/(a+b)
	km_inf = addnoise(a*tau_km)

	m_exp = 1 - exp(-dt/tau_m)
	h_exp = 1 - exp(-dt/tau_h)
	n_exp = 1 - exp(-dt/tau_n)
	km_exp= 1 - exp(-dt/tau_km)	

	
   UNITSON
}

FUNCTION vtrap(x,y) {  
    if (fabs(exp(x/y) - 1) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}
FUNCTION addnoise(input){
	addnoise=input
	if(NF>0){
		addnoise=input*normrand(1,NF*input*(1-input))
	}
	if(addnoise<0){addnoise=0}
	if(addnoise>1){addnoise=1}
	

}