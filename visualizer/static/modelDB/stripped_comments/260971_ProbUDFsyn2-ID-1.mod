NEURON {

        POINT_PROCESS ProbUDFsyn2  
        RANGE tau_r, tau_d
        RANGE Use, u, Dep, Fac, u0
        RANGE i, g, e, gmax
        NONSPECIFIC_CURRENT i
	POINTER rng
}

PARAMETER {

        tau_r = 0.2   (ms)  
        tau_d = 1.7    (ms)  
        Use = 1.0   (1)   
        Dep = 100   (ms)  
        Fac = 10   (ms)  
        e = 0     (mV)  
    	gmax = .001 (uS) 
    	u0 = 0 
}


   
VERBATIM

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);

ENDVERBATIM
  

ASSIGNED {

        v (mV)
        i (nA)
	g (uS)
        factor
	rng
	weight_NMDA
}

STATE {
        A       
        B       
}

INITIAL{

  LOCAL tp
        
	A = 0
  B = 0
	
        
	tp = (tau_r*tau_d)/(tau_d-tau_r)*log(tau_d/tau_r) 
	      
	factor = -exp(-tp/tau_r)+exp(-tp/tau_d) 
        factor = 1/factor
 
}

BREAKPOINT {

        SOLVE state METHOD cnexp
        g = gmax*(B-A) 
        i = g*(v-e) 
}

DERIVATIVE state{

        A' = -A/tau_r
        B' = -B/tau_d
}


NET_RECEIVE (weight, Pv, Pv_tmp, Pr, u, tsyn (ms)){
	
        INITIAL{
                Pv=1
                u=u0
                tsyn=t
            }

        
        if (Fac > 0) {
                u = u*exp(-(t - tsyn)/Fac) 
           } else {
                  u = Use  
           } 
           if(Fac > 0){
                  u = u + Use*(1-u) 
           }    

        
            Pv_tmp  = 1 - (1-Pv) * exp(-(t-tsyn)/Dep) 
                                                      
            Pr  = u * Pv_tmp                          
            Pv_tmp  = Pv_tmp - u * Pv_tmp             
            
            
                
		   if (erand() < Pr){
                    tsyn = t
	            Pv = Pv_tmp
                    A = A + weight*factor
                    B = B + weight*factor
                }
}

PROCEDURE setRNG() {
VERBATIM
    {
        
        void** pv = (void**)(&_p_rng);
        if( ifarg(1)) {
            *pv = nrn_random_arg(1);
        } else {
            *pv = (void*)0;
        }
    }
ENDVERBATIM
}

FUNCTION erand() {
VERBATIM
	    
        double value;
        if (_p_rng) {
                
                value = nrn_random_pick(_p_rng);
		        
                
                
                return value;
        }else{
ENDVERBATIM
                
                
                
                
                erand = exprand(1)
VERBATIM
        }
ENDVERBATIM
        
                       
                       
}