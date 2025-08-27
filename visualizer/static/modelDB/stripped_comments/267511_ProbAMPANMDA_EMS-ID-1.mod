NEURON {
    THREADSAFE
    POINT_PROCESS ProbAMPANMDA_EMS

    GLOBAL tau_r_AMPA
    RANGE tau_d_AMPA, g_AMPA, i_AMPA

    GLOBAL tau_r_NMDA, tau_d_NMDA
    RANGE g_NMDA, i_NMDA

    RANGE Use, u, Dep, Fac, u0, mg, tsyn
    RANGE unoccupied, occupied, Nrrp

    RANGE g, NMDA_ratio
    RANGE A_AMPA_step, B_AMPA_step, A_NMDA_step, B_NMDA_step
    GLOBAL e
    NONSPECIFIC_CURRENT i
    BBCOREPOINTER rng
    RANGE synapseID, selected_for_report, verboseLevel
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
        gmax = .001 (uS) 
        u0 = 0 
        Nrrp = 1 (1)  
        synapseID = 0
        verboseLevel = 0
        selected_for_report = 0
        NMDA_ratio = 0.71 (1) 
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
        mggate
        usingR123            

        
        unoccupied (1) 
        occupied   (1) 
        tsyn (ms) 
        u (1) 
}


STATE {

        A_AMPA       
        B_AMPA       
        A_NMDA       
        B_NMDA       
}


INITIAL {
        LOCAL tp_AMPA, tp_NMDA

        tsyn = 0
        u=u0

        
        unoccupied = 0
        occupied = Nrrp

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

        VERBATIM
        if( usingR123 ) {
            nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
        }
        ENDVERBATIM
}

BREAKPOINT {

        SOLVE state
        mggate = 1 / (1 + exp(0.062 (/mV) * -(v)) * (mg / 3.57 (mM))) 
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


NET_RECEIVE (weight,weight_AMPA, weight_NMDA, Psurv) {
    LOCAL result, ves, occu
    weight_AMPA = weight
    weight_NMDA = weight * NMDA_ratio
    
    

    INITIAL {
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
        A_AMPA = A_AMPA + ves/Nrrp*weight_AMPA*factor_AMPA
        B_AMPA = B_AMPA + ves/Nrrp*weight_AMPA*factor_AMPA
        A_NMDA = A_NMDA + ves/Nrrp*weight_NMDA*factor_NMDA
        B_NMDA = B_NMDA + ves/Nrrp*weight_NMDA*factor_NMDA

        if ( verboseLevel > 0 ) {
            UNITSOFF
            printf( "Release! %f at time %g
            UNITSON
        }

    } else {
        
        if ( verboseLevel > 0 ) {
            UNITSOFF
            printf( " || SYN_ID
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
    verboseLevel = 1-verboseLevel
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