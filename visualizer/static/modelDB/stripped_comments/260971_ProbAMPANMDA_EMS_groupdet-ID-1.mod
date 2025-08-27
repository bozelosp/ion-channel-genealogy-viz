NEURON {
    THREADSAFE
        POINT_PROCESS ProbAMPANMDA_EMS_groupdet
        RANGE tau_r_AMPA, tau_d_AMPA, tau_r_NMDA, tau_d_NMDA, Nsyns, Nevents, eventCounter
        RANGE Use, u, Dep, Fac, u0, mg, Rstate, tsyn_fac, u
        RANGE i, i_AMPA, i_NMDA, g_AMPA, g_NMDA, g, e, NMDA_ratio
        RANGE A_AMPA_step, B_AMPA_step, A_NMDA_step, B_NMDA_step
        NONSPECIFIC_CURRENT i
        POINTER rng
        RANGE synapseID, verboseLevel
}

PARAMETER {


        tau_r_AMPA = 0.2   (ms)  
        tau_d_AMPA = 1.7    (ms)  
        tau_r_NMDA = 0.29   (ms) 
        tau_d_NMDA = 43     (ms) 
        Use = 1.0   (1)   
        Dep = 100   (ms)  
        Fac = 10   (ms)  
        e = 0     (mV)  
        mg = 1   (mM)  
        mggate
        gmax = .001 (uS) 
        u0 = 0 
        synapseID = 0
        verboseLevel = 0
	NMDA_ratio = 0.71 (1) 
	Nsyns = 10
        Nevents = 0 
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
        i_AMPA (nA)
        i_NMDA (nA)
        g_AMPA (uS)
        g_NMDA (uS)
        g (uS)
        factor_AMPA
        factor_NMDA
        A_AMPA_step
        B_AMPA_step
        A_NMDA_step
        B_NMDA_step
        rng

	
	
        
        
        
	Rstate (1) 
	tsyn_fac (ms) 
	u (1) 
        space       
        space2       
        eventCounter 
}

STATE {

        A_AMPA       
        B_AMPA       
        A_NMDA       
        B_NMDA       
}

INITIAL{

        LOCAL tp_AMPA, tp_NMDA

	Rstate=1
	tsyn_fac=0
	u=u0
        
        A_AMPA = 0
        B_AMPA = 0
        
        A_NMDA = 0
        B_NMDA = 0
        
        tp_AMPA = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA) 
        tp_NMDA = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA) 
        
        factor_AMPA = -exp(-tp_AMPA/tau_r_AMPA)+exp(-tp_AMPA/tau_d_AMPA) 
        factor_AMPA = 1/factor_AMPA
        
        factor_NMDA = -exp(-tp_NMDA/tau_r_NMDA)+exp(-tp_NMDA/tau_d_NMDA) 
        factor_NMDA = 1/factor_NMDA

        A_AMPA_step = exp(dt*(( - 1.0 ) / tau_r_AMPA))
        B_AMPA_step = exp(dt*(( - 1.0 ) / tau_d_AMPA))
        A_NMDA_step = exp(dt*(( - 1.0 ) / tau_r_NMDA))
        B_NMDA_step = exp(dt*(( - 1.0 ) / tau_d_NMDA))
        eventCounter = 0
}

BREAKPOINT {

        SOLVE state
        mggate = 1 / (1 + (mg/4.1 (mM))*exp(0.063 (/mV)*(-v)))
        g_AMPA = gmax*(B_AMPA-A_AMPA) 
        g_NMDA = gmax*(B_NMDA-A_NMDA) * mggate 
        g = g_AMPA + g_NMDA
        i_AMPA = g_AMPA*(v-e) 
        i_NMDA = g_NMDA*(v-e) 
        i = i_AMPA + i_NMDA
}

PROCEDURE state() {
        A_AMPA = A_AMPA*A_AMPA_step
        B_AMPA = B_AMPA*B_AMPA_step
        A_NMDA = A_NMDA*A_NMDA_step
        B_NMDA = B_NMDA*B_NMDA_step
}


NET_RECEIVE (weight, Psurv, myInd, tsyn (ms), tsyn_fac (ms), u){
        LOCAL result
	
	
	
	
	
	
        INITIAL{
                tsyn=t
                eventCounter=0
        }

        
        VERBATIM
          void** vv = (void**)(&space);
          void** vv2 = (void**)(&space2);
          double *x;
          int nx = vector_instance_px(*vv, &x);
          double *x2;
          int nx2 = vector_instance_px(*vv2, &x2);
          int myInd = 0;
          if (eventCounter < nx2) {
            myInd = x2[(int)eventCounter];
            
          }
          else printf("eventCounter >= nx2! t = %lf, eventCounter = %lf, nx2 = %i\n",t, eventCounter, nx2);
          eventCounter++;
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
                         A_AMPA = A_AMPA + weight*factor_AMPA
                         B_AMPA = B_AMPA + weight*factor_AMPA
                         A_NMDA = A_NMDA + weight*NMDA_ratio*factor_NMDA
                         B_NMDA = B_NMDA + weight*NMDA_ratio*factor_NMDA
                         
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
    verboseLevel = 1-verboseLevel
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

PROCEDURE setVec2() {    
  VERBATIM
  void** vv;
  vv = (void**)(&space2);
  *vv = (void*)0;
  if (ifarg(1)) {
    *vv = vector_arg(1);
    Nevents = vector_capacity(*vv);
  }
  ENDVERBATIM
}