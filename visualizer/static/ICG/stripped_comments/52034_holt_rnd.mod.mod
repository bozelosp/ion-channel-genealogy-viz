NEURON {
        SUFFIX random
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

VERBATIM
double holt_random();
void holt_seed();
double holt_normrand();
double holt_exprand();
int holt_poisrand();
ENDVERBATIM




FUNCTION Uniform() {
  Uniform = holt_random()	
}





FUNCTION Poisson(mean) {
VERBATIM
  _lPoisson = holt_poisrand(_lmean);
ENDVERBATIM
}




FUNCTION Exp() {
  Exp = holt_exprand()		
}




FUNCTION Normal(mean, std_dev) {
VERBATIM
  _lNormal = holt_normrand(_lmean, _lstd_dev); 
ENDVERBATIM
}







PROCEDURE Seed(seedval, algorithm) {
VERBATIM
  holt_seed((int)_lseedval, (int)_lalgorithm);
ENDVERBATIM
}

VERBATIM


int seed = 0;			
#define DRAND48 0		
#ifndef hpux			
#define RANDOM 1
#define N_GENERATORS 2		
#else
#define N_GENERATORS 1
#endif
#define DEFAULT_GENERATOR DRAND48 
int random_generator;		

#ifdef hpux
extern double drand48();
#endif
#ifndef __alpha
extern long random();
extern double drand48();
#endif



void
holt_seed(seedval, gen_code)
int seedval, gen_code;
{
  if (seedval == 0)		
	seedval=3491;


  seed = seedval;

  if (gen_code >= N_GENERATORS)	
    gen_code = DEFAULT_GENERATOR; 
  random_generator = gen_code;	

  switch (random_generator)	
  {				
  case DRAND48:    srand48((long)seedval);    break;



  }
}


double
holt_random()
{
  if (seed == 0)		
    holt_seed(0, DEFAULT_GENERATOR); 

  switch (random_generator)	
  {
  case DRAND48:   return drand48();



  }
}
























double 
holt_normrand(mean, std_dev)
double mean, std_dev;
{
    double s, v1, v2;
    extern double sqrt(), log();

    s = 1.0;
    while (s >= 1.0)
    {
	v1 = 2.0 * holt_random() - 1.0;
	v2 = 2.0 * holt_random() - 1.0;
	s = (v1 * v1) + (v2 * v2);
    }
    v2 = v1 * sqrt(-2.0 * log(s) / s);
    return (v2 * std_dev + mean);
}


#define PI 3.141592654

int
holt_poisrand(xm)
double xm;
{
  float gammln();
  double drand48();
  static float sq,alxm,g,oldm=(-1.0);
  double em,tt,y;
  double exp();

  if (xm < 12.0) {
    if (xm != oldm) {
      oldm=xm;
      g=exp(-xm);
    }
    em = -1;
    tt=1.0;
    do {
      ++em;
      
      tt *= holt_random();
    } while (tt > g);
  } else {
    if (xm != oldm) {
      oldm=xm;
      sq=sqrt(2.0*xm);
      alxm=log(xm);
      g=xm*alxm-gammln(xm+1.0);
    }
    do {
      do {
	
	y=tan(PI*holt_random());
	em=sq*y+xm;
      } while (em < 0.0);
      em=floor(em);
      tt=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
    } while (holt_random() > tt);
  }
  return (int)em;
}
#undef PI


float
gammln(xx)
float xx;
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























double 
holt_exprand()
{
    extern double log();

    return (-log(holt_random()));
}
ENDVERBATIM