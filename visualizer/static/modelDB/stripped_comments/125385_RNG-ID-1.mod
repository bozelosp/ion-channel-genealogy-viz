INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
    init_seed = -1           
                        
                        
    DONT_VECTORIZE          
}

NEURON {
    POINT_PROCESS RNG
    RANGE init_seed
    GLOBAL DONT_VECTORIZE   
}

INITIAL {
VERBATIM
    
    if (init_seed != -1)
        set_seed(init_seed);
    
ENDVERBATIM
}

VERBATIM

#include <math.h>
#include <limits.h>

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



static long idum_RNG = -1;
static long iy_RNG = 0; 
static long iv_RNG[r_ntab];



static void
set_seed(long seed) {
    int j, k;
    
    idum_RNG = fabs(seed);
    fprintf(stderr,"RNG: seed is set to %ld\n",idum_RNG); 

    for ( j = r_ntab + 7; j>=0;j--){
        k = (idum_RNG/r_iq);
        idum_RNG = r_ia*(idum_RNG-k*r_iq)-r_ir*k;
        if (idum_RNG<0) idum_RNG += r_im;
        if (j<r_ntab) iv_RNG[j] = idum_RNG;
    }
    iy_RNG = iv_RNG[0];
}




static double
uniran1()
{     
    int  j;
     long k;
     double temp; 
    k = (idum_RNG/r_iq);
     idum_RNG = r_ia*(idum_RNG-k*r_iq)-r_ir*k;
     if (idum_RNG<0) idum_RNG += r_im;
     j = iy_RNG/r_ndiv;
     iy_RNG = iv_RNG[j];
     iv_RNG[j] = idum_RNG;
     if ((temp = r_am*iy_RNG) >r_rnmx) return r_rnmx;
     else return temp;
}



static double
gasdev() {
     double fac,r,v1,v2;
    static int iset_RNG = 0;
    static double gset_RNG;

     if  (iset_RNG == 0) {
        do {
            v1=2.0*uniran1()-1.0;
               v2=2.0*uniran1()-1.0;
               r=v1*v1+v2*v2;
        } while (r >= 1.0 || r == 0.0);
        fac=sqrt(-2.0*log(r)/r);
        gset_RNG=v1*fac;
        iset_RNG=1;
        return v2*fac;
    } else {
        iset_RNG=0;
          return gset_RNG;
    }
}



double
bnldev(double ppr, int nnr) {
    int j;
    static int nold=(-1);
    double am,em,g,angle,p,bnl,sq,bt,y;
    static double pold=(-1.0),pc,plog,pclog,en,oldg;
    
    
     
    
    p=(ppr <= 0.5 ? ppr : 1.0-ppr);
    am=nnr*p;
    if (nnr < 25) {
        bnl=0.0;
        for (j=1;j<=nnr;j++)
            if (uniran1() < p) bnl += 1.0;
    }
    else if (am < 1.0) {
        g=exp(-am);
        bt=1.0;
        for (j=0;j<=nnr;j++) {
            bt *= uniran1();
            if (bt < g) break;
        }
        bnl=(j <= nnr ? j : nnr);
    }
    else {
        if (nnr != nold) {
            en=nnr;
            oldg=gammln(en+1.0);
            nold=nnr;
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
                angle=PI*uniran1();
                    angle=PI*uniran1();
                y=tan(angle);
                em=sq*y+am;
            } while (em < 0.0 || em >= (en+1.0));
            em=floor(em);
                bt=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0) - 
                gammln(en-em+1.0)+em*plog+(en-em)*pclog);
        } while (uniran1() > bt);
        bnl=em;
    }
    if (p != ppr) bnl=nnr-bnl;
    
    
   
    
    return bnl;
}

#undef PI 
#undef r_ia 
#undef r_im 
#undef r_am 
#undef r_iq  
#undef r_ir 
#undef r_ntab 
#undef r_ndiv 
#undef r_eps 
#undef r_rnmx 
ENDVERBATIM



PROCEDURE SetSeedNow(seed) {
VERBATIM
    set_seed((long)_lseed); 
ENDVERBATIM
}



FUNCTION UniDev () {
VERBATIM
    _lUniDev = uniran1();   
ENDVERBATIM
}
    


FUNCTION BnlDev (ppr, nnr) {
VERBATIM
    _lBnlDev = bnldev(_lppr, (int)_lnnr);  
ENDVERBATIM
}



PROCEDURE test_ran1() {
    LOCAL j,x
      VERBATIM 
      FILE *out; 
      out = fopen("test_ran1.out","w");
      ENDVERBATIM 

    j = 0
    while (j<10000) {
      VERBATIM 
       fprintf(out,"%lg\n",uniran1());
      ENDVERBATIM 
      j = j + 1
    }
      VERBATIM 
      fclose(out);
      ENDVERBATIM 
   
}       



PROCEDURE test_gasdev() {
    LOCAL j,x
      VERBATIM 
      FILE *out; 
      out = fopen("test_gasdev.out","w");
      ENDVERBATIM 

    j = 0
    while (j<10000) {
      VERBATIM 
       fprintf(out,"%lg\n",gasdev());
      ENDVERBATIM 
      j = j + 1
    }
      VERBATIM 
      fclose(out);
      ENDVERBATIM 
   
}