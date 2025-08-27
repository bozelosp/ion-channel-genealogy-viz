INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX StochKv
    THREADSAFE
    USEION k READ ek WRITE ik
    RANGE N,eta, gk, gamma, deterministic, gkbar, ik
    RANGE ninf, ntau,a,b,P_a,P_b
    GLOBAL Ra, Rb
    GLOBAL vmin, vmax, q10, temp
    RANGE tadj
    BBCOREPOINTER rng
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (pS) = (picosiemens)
    (S) = (siemens)
    (um) = (micron)
} 

PARAMETER {
    v           (mV)
    dt      (ms)
    area    (um2)
    
    gamma  =  30          (pS)
    eta              (1/um2)
    gkbar = .75      (S/cm2)
    
    tha  = -40   (mV)        
    qa   = 9            
    Ra   = 0.02 (/ms)       
    Rb   = 0.002    (/ms)       
    
    celsius (degC)
    temp = 23 (degC)   
    q10 = 2.3               
    
    deterministic = 0   
    vmin = -120 (mV)    
    vmax = 100  (mV)
} 

ASSIGNED {
    a       (/ms)
    b       (/ms)
    ik      (mA/cm2)
    gk      (S/cm2)
    ek      (mV)
    ninf        
    ntau (ms)   
    tadj

    N 
    scale_dens (pS/um2) 
    P_a     
    P_b     

    rng
    usingR123

    n0_n1_new

}


STATE {
    n         
    N0 N1       
    n0_n1 n1_n0 
}


   
VERBATIM
#include "nrnran123.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef CORENEURON_BUILD
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#endif

ENDVERBATIM


INITIAL {
    VERBATIM
    if( usingR123 ) {
        nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
    }
    ENDVERBATIM
  
    eta = gkbar / gamma
    trates(v)
    n = ninf
    scale_dens = gamma/area
    N = floor(eta*area + 0.5)
    
    N1 = floor(n * N + 0.5)
    N0 = N-N1       
    
    n0_n1 = 0
    n1_n0 = 0
}



BREAKPOINT {
  SOLVE states
  
  gk =  (strap(N1) * scale_dens * tadj)
  
  ik = 1e-4 * gk * (v - ek)
} 




PROCEDURE states() {

    trates(v)
    
    P_a = strap(a*dt)
    P_b = strap(b*dt)

    
    ChkProb( P_a)
    ChkProb( P_b)
    
    
    n0_n1 = BnlDev(P_a, N0)
    n1_n0 = BnlDev(P_b, N1)

    
    N0    = strap(N0 - n0_n1 + n1_n0)
    N1    = N - N0
}



PROCEDURE trates(v (mV)) {     
    TABLE ntau, ninf, a, b, tadj
    DEPEND dt, Ra, Rb, tha, qa, q10, temp, celsius
    FROM vmin TO vmax WITH 199
    
    tadj = q10 ^ ((celsius - temp)/(10 (K)))
    a = SigmoidRate(v, tha, Ra, qa)
    a = a * tadj
    b = SigmoidRate(-v, -tha, Rb, qa)
    b = b * tadj
    ntau = 1/(a+b)
    ninf = a*ntau
}





FUNCTION SigmoidRate(v (mV),th (mV),a (1/ms),q) (1/ms){
    UNITSOFF
    if (fabs(v-th) > 1e-6 ) {
        SigmoidRate = a * (v - th) / (1 - exp(-(v - th)/q))
    UNITSON

    } else {
        SigmoidRate = a * q
    }
}   




FUNCTION strap(x) {
    if (x < 0) {
        strap = 0
VERBATIM
        fprintf (stderr,"skv.mod:strap: negative state");
ENDVERBATIM
    } else {
        strap = x
    }
}



PROCEDURE ChkProb(p) {

  if (p < 0.0 || p > 1.0) {
    VERBATIM


    ENDVERBATIM
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
    double value;
    if( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
#ifndef CORENEURON_BUILD
        value = nrn_random_pick(_p_rng);
#endif
    } else {
        value = 0.5;
    }
    _lurand = value;
ENDVERBATIM
}


FUNCTION brand(P, N) {

VERBATIM
        

        
        double value = 0.0;
        int i;
        for (i = 0; i < _lN; i++) {
           if (urand(_threadargs_) < _lP) {
              value = value + 1;
           }
        }
        return(value);

ENDVERBATIM

        brand = value
}

VERBATIM
#define        PI 3.141592654
#define        r_ia     16807
#define        r_im     2147483647
#define        r_am     (1.0/r_im)
#define        r_iq     127773
#define        r_ir     2836
#define        r_ntab   32
#define        r_ndiv   (1+(r_im-1)/r_ntab)
#define        r_eps    1.2e-7
#define        r_rnmx   (1.0-r_eps)
ENDVERBATIM

VERBATIM


static double
gammln(double xx)
{
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
        -1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    x=xx-1.0;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.0;
    for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
    }
    return -tmp+log(2.50662827465*ser);
}
ENDVERBATIM




FUNCTION BnlDev (ppr, nnr) {

VERBATIM
        int j;
        static int nold=(-1);
        double am,em,g,angle,p,bnl,sq,bt,y;
        static double pold=(-1.0),pc,plog,pclog,en,oldg;
        
        
         
        
        p=(_lppr <= 0.5 ? _lppr : 1.0-_lppr);
        am=_lnnr*p;
        if (_lnnr < 25) {
            bnl=0.0;
            for (j=1;j<=_lnnr;j++)
                if (urand(_threadargs_) < p) bnl += 1.0;
        }
        else if (am < 1.0) {
            g=exp(-am);
            bt=1.0;
            for (j=0;j<=_lnnr;j++) {
                bt *= urand(_threadargs_);
                if (bt < g) break;
            }
            bnl=(j <= _lnnr ? j : _lnnr);
        }
        else {
            if (_lnnr != nold) {
                en=_lnnr;
                oldg=gammln(en+1.0);
                nold=_lnnr;
            }
            if (p != pold) {
                pc=1.0-p;
                 plog=log(p);
                pclog=log(pc);
                pold=p;
            }
            sq=sqrt(2.0*am*pc);
            do {
                do {
                    angle=PI*urand(_threadargs_);
                        angle=PI*urand(_threadargs_);
                    y=tan(angle);
                    em=sq*y+am;
                } while (em < 0.0 || em >= (en+1.0));
                em=floor(em);
                    bt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) - 
                    gammln(en-em+1.0)+em*plog+(en-em)*pclog);
            } while (urand(_threadargs_) > bt);
            bnl=em;
        }
        if (p != _lppr) bnl=_lnnr-bnl;
        
        
       
        
        return bnl;
        
    ENDVERBATIM
    BnlDev = bnl
}  

VERBATIM
static void bbcore_write(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    if (d) {
        uint32_t* di = ((uint32_t*)d) + *offset;
      
      if (!_p_rng) {
        di[0] = 0; di[1] = 0, di[2] = 0;
      }else{
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        nrnran123_getids3(*pv, di, di+1, di+2);
        
        char which;
        nrnran123_getseq(*pv, di+3, &which);
        di[4] = (int)which;
      }
      
    }
    *offset += 5;
}
static void bbcore_read(double* x, int* d, int* xx, int* offset, _threadargsproto_) {
    assert(!_p_rng);
    uint32_t* di = ((uint32_t*)d) + *offset;
        if (di[0] != 0 || di[1] != 0|| di[2] != 0)
        {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
      *pv = nrnran123_newstream3(di[0], di[1], di[2]);
      nrnran123_setseq(*pv, di[3], (char)di[4]);
        }
      
    *offset += 5;
}
ENDVERBATIM

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
                
                if (*xdir == -1.) {
                    if( usingR123 ) {
                        *xdir = 2.0;
                    } else {
                        *xdir = 1.0;
                    }
                    return 0.0;
                }
                else if (*xdir == 0.) {
                    if( usingR123 ) {
                        uint32_t seq;
                        char which;
                        nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );
                        xval[0] = (double) seq;
                        xval[1] = (double) which;
                    } else {
                        xval[0] = (double)nrn_get_random_sequence(_p_rng);
                    }
                } else{
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