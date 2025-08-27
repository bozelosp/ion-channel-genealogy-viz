NEURON {

        POINT_PROCESS NMDA
        RANGE  tau_r_NMDA, tau_d_NMDA,n_NMDA,gama_NMDA
        RANGE Use
        RANGE i,  i_NMDA,  g_NMDA, e, gmax
        NONSPECIFIC_CURRENT i
}

PARAMETER {

    	n_NMDA = 0.28011 (/mM)	
    	gama_NMDA = 0.062 (/mV) 
	   tau_r_NMDA = 0.3   (ms) 
        tau_d_NMDA = 43     (ms) 
        Use = 1.0   (1)   

        e = 0     (mV)  
	    mg = 1   (mM)  
        mggate
    	
    	u0 = 0 
}


   

  

ASSIGNED {

        v (mV)
        i (nA)
	i_NMDA (nA)
	g_NMDA (uS)
	factor_NMDA
	
}

STATE {

       
	A_NMDA       
    B_NMDA       
}

INITIAL{

    LOCAL  tp_NMDA
	
	A_NMDA = 0
	B_NMDA = 0
        
	tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA) 
        
	
	
	factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA) 
    factor_NMDA = 1/factor_NMDA
   
}

BREAKPOINT {

    SOLVE state METHOD cnexp
	mggate = 1 / (1 + exp(gama_NMDA  * -(v)) * (n_NMDA)) 
	g_NMDA = (B_NMDA-A_NMDA) * mggate 
       
	i_NMDA = g_NMDA*(v-e) 
	i =  i_NMDA
}

DERIVATIVE state{

    
	A_NMDA' = -A_NMDA/tau_r_NMDA
    B_NMDA' = -B_NMDA/tau_d_NMDA
}





NET_RECEIVE (weight, weight_NMDA){
	
	weight_NMDA = weight
	

     
	A_NMDA = A_NMDA + weight_NMDA*factor_NMDA
    B_NMDA = B_NMDA + weight_NMDA*factor_NMDA

              
}