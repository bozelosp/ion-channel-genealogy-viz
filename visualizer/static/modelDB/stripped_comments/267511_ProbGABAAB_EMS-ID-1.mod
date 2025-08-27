NEURON {
    THREADSAFE
	POINT_PROCESS ProbGABAAB_EMS
	RANGE tau_r_GABAA, tau_d_GABAA, tau_r_GABAB, tau_d_GABAB
	RANGE Use, u, Dep, Fac, u0, tsyn
    RANGE unoccupied, occupied, Nrrp
	RANGE i,i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, e_GABAA, e_GABAB, GABAB_ratio
        RANGE A_GABAA_step, B_GABAA_step, A_GABAB_step, B_GABAB_step
	NONSPECIFIC_CURRENT i
    BBCOREPOINTER rng
    RANGE synapseID, selected_for_report, verboseLevel
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
    Nrrp = 1 (1)  
    synapseID = 0
    verboseLevel = 0
    selected_for_report = 0
	GABAB_ratio = 0 (1) 
}



VERBATIM
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include "nrnran123.h"

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
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
        usingR123            

    
    unoccupied (1) 
    occupied   (1) 
    tsyn (ms) 
    u (1) 
}

STATE {
        A_GABAA       
        B_GABAA       
        A_GABAB       
        B_GABAB       
}

INITIAL{
        LOCAL tp_GABAA, tp_GABAB

        tsyn = 0
        u=u0

        
        unoccupied = 0
        occupied = Nrrp

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

        VERBATIM
        if( usingR123 ) {
            nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
        }
        ENDVERBATIM
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


NET_RECEIVE (weight, weight_GABAA, weight_GABAB, Psurv){
    LOCAL result, ves, occu
    weight_GABAA = weight
    weight_GABAB = weight*GABAB_ratio
    
    


    INITIAL{
    }

    
    
    
    if(  weight <= 0 || t < 0 ) {
VERBATIM
        return;
ENDVERBATIM
    }

    
    if (Fac > 0) {
            u = u*exp(-(t - tsyn)/Fac) 
       } else {
              u = Use
       }
       if(Fac > 0){
              u = u + Use*(1-u) 
       }

    
    FROM counter = 0 TO (unoccupied - 1) {
        
        Psurv = exp(-(t-tsyn)/Dep)
        result = urand()
        if (result>Psurv) {
            occupied = occupied + 1     
            if( verboseLevel > 0 ) {
                UNITSOFF
                printf( "Recovered! %f at time %g
                UNITSON
            }
        }
    }

    ves = 0                  
    occu = occupied - 1  

    FROM counter = 0 TO occu {
        
        result = urand()
        if (result<u) {
            
            occupied = occupied - 1  
            ves = ves + 1            
        }
    }

    
    unoccupied = Nrrp - occupied

    
    
    
    
    tsyn = t

    if (ves > 0) { 
        A_GABAA = A_GABAA + ves/Nrrp*weight_GABAA*factor_GABAA
        B_GABAA = B_GABAA + ves/Nrrp*weight_GABAA*factor_GABAA
        A_GABAB = A_GABAB + ves/Nrrp*weight_GABAB*factor_GABAB
        B_GABAB = B_GABAB + ves/Nrrp*weight_GABAB*factor_GABAB

        if( verboseLevel > 0 ) {
            UNITSOFF
            printf( "Release! %f at time %g
            UNITSON
        }

    } else {
        
        if ( verboseLevel > 0 ) {
            UNITSOFF
            printf("Failure! %f at time %g
            UNITSON
        }
    }

}


PROCEDURE setRNG() {
VERBATIM
    #ifndef CORENEURON_BUILD
    
    
    
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        uint32_t a2 = 0;
        uint32_t a3 = 0;

        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            a3 = (uint32_t)*getarg(3);
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);
        usingR123 = 1;
    } else if( ifarg(1) ) {   
        void** pv = (void**)(&_p_rng);
        *pv = nrn_random_arg(1);
    } else {  
        void** pv = (void**)(&_p_rng);
        *pv = (void*)0;
    }
    #endif
ENDVERBATIM
}


FUNCTION urand() {
VERBATIM
    double value = 0.0;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
        #ifndef CORENEURON_BUILD
        value = nrn_random_pick(_p_rng);
        #endif
    } else {
        
        value = 0.0;
    }
    _lurand = value;
ENDVERBATIM
}


FUNCTION bbsavestate() {
        bbsavestate = 0
VERBATIM
#ifndef CORENEURON_BUILD
        
        
        double *xdir, *xval, *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng) {
            
            if (*xdir == -1) {  
                if( usingR123 ) {
                    *xdir = 2.0;
                } else {
                    *xdir = 1.0;
                }
                return 0.0;
            } else if(*xdir ==0 ) {  
                if( usingR123 ) {
                    uint32_t seq;
                    char which;
                    nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );
                    xval[0] = (double) seq;
                    xval[1] = (double) which;
                } else {
                    xval[0] = (double)nrn_get_random_sequence(_p_rng);
                }
            } else {  
                if( usingR123 ) {
                    nrnran123_setseq( (nrnran123_State*)_p_rng, (uint32_t)xval[0], (char)xval[1] );
                } else {
                    nrn_set_random_sequence(_p_rng, (long)(xval[0]));
                }
            }
        }
#endif
ENDVERBATIM
}


FUNCTION toggleVerbose() {
    verboseLevel = 1 - verboseLevel
}


VERBATIM
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
   if (d) {
    
    uint32_t* di = ((uint32_t*)d) + *offset;
    nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
    nrnran123_getids3(*pv, di, di+1, di+2);

    
    char which;
    nrnran123_getseq(*pv, di+3, &which);
    di[4] = (int)which;
    
   }
  *offset += 5;
}

static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
  assert(!_p_rng);
  uint32_t* di = ((uint32_t*)d) + *offset;
  if (di[0] != 0 || di[1] != 0 || di[2] != 0) {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
      *pv = nrnran123_newstream3(di[0], di[1], di[2]);

      
      unsigned char which = (unsigned char)di[4];
      nrnran123_setseq(*pv, di[3], which);
  }
  
  *offset += 5;
}
ENDVERBATIM