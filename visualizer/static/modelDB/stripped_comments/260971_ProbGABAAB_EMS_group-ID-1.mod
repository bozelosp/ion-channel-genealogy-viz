NEURON {
    THREADSAFE
	POINT_PROCESS ProbGABAAB_EMS_group
	RANGE tau_r_GABAA, tau_d_GABAA, tau_r_GABAB, tau_d_GABAB, Nsyns
	RANGE Use, u, Dep, Fac, u0, Rstate, tsyn_fac, u
	RANGE i,i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, e_GABAA, e_GABAB, GABAB_ratio, gmax
        RANGE A_GABAA_step, B_GABAA_step, A_GABAB_step, B_GABAB_step
	NONSPECIFIC_CURRENT i
    POINTER rng
    RANGE synapseID, verboseLevel
}

PARAMETER {
	tau_r_GABAA  = 0.2   (ms)  
	tau_d_GABAA = 8   (ms)  
    tau_r_GABAB  = 3.5   (ms)  
	tau_d_GABAB = 260.9   (ms)  
	Use        = 1.0   (1)   
	Dep   = 100   (ms)  
	Fac   = 10   (ms)  
	e_GABAA    = -80     (mV)  
    e_GABAB    = -97     (mV)  
    gmax = .001 (uS) 
    u0 = 0 
    synapseID = 0
    verboseLevel = 0
	GABAB_ratio = 0 (1) 
        Nsyns = 10
}


   
VERBATIM
#include<stdlib.h>
#include<stdio.h>
#include<math.h>

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);

extern int ifarg(int iarg);
extern int vector_capacity(void* vv);
extern void* vector_arg(int iarg);
ENDVERBATIM
  

ASSIGNED {
	v (mV)
	i (nA)
        i_GABAA (nA)
        i_GABAB (nA)
        g_GABAA (uS)
        g_GABAB (uS)
        A_GABAA_step
        B_GABAA_step
        A_GABAB_step
        B_GABAB_step
	g (uS)
	factor_GABAA
        factor_GABAB
        rng

       
       
	 
	 
	 
       Rstate (1) 
       tsyn_fac (ms) 
       u (1) 
       space       

}

STATE {
        A_GABAA       
        B_GABAA       
        A_GABAB       
        B_GABAB       
}

INITIAL{

        LOCAL tp_GABAA, tp_GABAB

	Rstate=1
	tsyn_fac=0
	u=u0
        
        A_GABAA = 0
        B_GABAA = 0
        
        A_GABAB = 0
        B_GABAB = 0
        
        tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA) 
        tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB) 
        
        factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) 
        factor_GABAA = 1/factor_GABAA
        
        factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) 
        factor_GABAB = 1/factor_GABAB
        
        A_GABAA_step = exp(dt*(( - 1.0 ) / tau_r_GABAA))
        B_GABAA_step = exp(dt*(( - 1.0 ) / tau_d_GABAA))
        A_GABAB_step = exp(dt*(( - 1.0 ) / tau_r_GABAB))
        B_GABAB_step = exp(dt*(( - 1.0 ) / tau_d_GABAB))
}

BREAKPOINT {
	SOLVE state
	
        g_GABAA = gmax*(B_GABAA-A_GABAA) 
        g_GABAB = gmax*(B_GABAB-A_GABAB) 
        g = g_GABAA + g_GABAB
        i_GABAA = g_GABAA*(v-e_GABAA) 
        i_GABAB = g_GABAB*(v-e_GABAB) 
        i = i_GABAA + i_GABAB
}

PROCEDURE state() {
        A_GABAA = A_GABAA*A_GABAA_step
        B_GABAA = B_GABAA*B_GABAA_step
        A_GABAB = A_GABAB*A_GABAB_step
        B_GABAB = B_GABAB*B_GABAB_step
}


NET_RECEIVE (weight, Psurv, myInd, tsyn (ms), tsyn_fac (ms), u){
    LOCAL result
    
    
    
    
    
    INITIAL{
		tsyn=t
    }

    
    VERBATIM
      void** vv = (void**)(&space);
      double *x;
      int nx = vector_instance_px(*vv, &x);
      int myInd = rand()%((int)Nsyns);
      _args[2] = myInd;
      _args[3] = x[myInd];                
      _args[4] = x[myInd+(int)Nsyns];     
      _args[5] = x[myInd+2*((int)Nsyns)]; 
    ENDVERBATIM

    
    if(  !(weight > 0) ) {
VERBATIM
        return;
ENDVERBATIM
    }

        
        if (Fac > 0) {
                u = u*exp(-(t - tsyn_fac)/Fac) 
           } else {
                  u = Use  
           } 
           if(Fac > 0){
                  u = u + Use*(1-u) 
           }    

	   
	   
	   tsyn_fac = t

           
	   if (Rstate == 0) {
	   
	          Psurv = exp(-(t-tsyn)/Dep)
		  result = urand()
		  if (result>Psurv) {
		         Rstate = 1     

                         if( verboseLevel > 0 ) {
                             printf( "Recovered! %f at time %g
                         }

		  }
		  else {
		         
		         tsyn = t
                         if( verboseLevel > 0 ) {
                             printf( "Failed to recover! %f at time %g
                         }
		  }
           }	   
	   
	   if (Rstate == 1) {
   	          result = urand()
		  if (result<u) {
		  
   		         tsyn = t
			 Rstate = 0

                         A_GABAA = A_GABAA + weight*factor_GABAA
                         B_GABAA = B_GABAA + weight*factor_GABAA
                         A_GABAB = A_GABAB + weight*GABAB_ratio*factor_GABAB
                         B_GABAB = B_GABAB + weight*GABAB_ratio*factor_GABAB
                         
                         if( verboseLevel > 0 ) {
                             printf( "Release! %f at time %g
                         }
		  		  
		  }
		  else {
		         if( verboseLevel > 0 ) {
			     printf("Failure! %f at time %g
		         }

		  }

	   }
        VERBATIM
          x[myInd] = _args[3];                
          x[myInd+(int)Nsyns] = _args[4];     
          x[myInd+2*((int)Nsyns)] = _args[5]; 
        ENDVERBATIM

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

FUNCTION urand() {
VERBATIM
        double value;
        if (_p_rng) {
                
                value = nrn_random_pick(_p_rng);
                
                return value;
        }else{
ENDVERBATIM
                
                
                
                
                value = scop_random(1)
VERBATIM
        }
ENDVERBATIM
        urand = value
}

FUNCTION toggleVerbose() {
    verboseLevel = 1 - verboseLevel
}

PROCEDURE setVec() {    
                        
  VERBATIM
  void** vv;
  vv = (void**)(&space);
  *vv = (void*)0;
  if (ifarg(1)) {
    *vv = vector_arg(1);
    Nsyns = vector_capacity(*vv)/3;
  }
  ENDVERBATIM
}