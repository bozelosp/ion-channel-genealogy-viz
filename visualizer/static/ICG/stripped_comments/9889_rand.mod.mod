INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX nothing
}

VERBATIM
#include <stdlib.h>
#include <math.h>

#include <limits.h> 


#if !defined(MAXLONG)
#include <limits.h>
#define MAXLONG LONG_MAX
#endif

extern double my_drand48();
extern void my_srand48();
#undef drand48
#undef srand48
#define drand48 my_drand48
#define srand48 my_srand48

extern double drand48();
#define random()                drand48()*MAXLONG
#define initstate(c1,c2,c3)     srand48(c1)

static long state2[32] = {
	470594912, 650447616, 310934240, 695012864, 850358912,
61088076, 481306752, 786902080, 224042800, 805177664, 938284096,
145937936, 622867968, 160207584, 977329216, 716234240, 127727624,
415316352, 870137472, 18664444, 330872224, 93728752, 914779200,
736261248, 643647616, 755802688, 213052336, 410240448, 218974736,
109419280, 178026128, 689569664
};
ENDVERBATIM



FUNCTION fseed(seed) {
VERBATIM
    initstate((unsigned)_lseed,(char *)state2,32);
	_lfseed = _lseed;
ENDVERBATIM
}


FUNCTION n_rand() { 
VERBATIM
    _ln_rand = ((double)random()) / (((double)MAXLONG) + 1.);
ENDVERBATIM
}


FUNCTION fran(l, h) { 
VERBATIM
{
	int low, high;
    double num, imax, *getarg();
    
	low = (int)_ll;
	high = (int)_lh;
    imax = high-low+1; 
    _lfran = (double)(low + (int) (imax*n_rand()));  
}
ENDVERBATIM
}


FUNCTION u_rand() { 
VERBATIM
    _lu_rand = (((double)random()) / ((double)MAXLONG));
ENDVERBATIM
}
    

FUNCTION norm() { 
VERBATIM
{
    static int iset = 0;
    static float gset;
    float fac, r , v1, v2;
    double sqrt();

    if (iset == 0) {
        do {
	    	v1 = 2.0 * n_rand() - 1.0;
		    v2 = 2.0 * n_rand() - 1.0;
		    r = v1 * v1 + v2 * v2;
	    } while (r >= 1.0);

        fac = (float)sqrt(-2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;
        _lnorm = v2 * fac;

    } else {
        iset = 0;
        _lnorm = (double)gset;
    }
}
ENDVERBATIM
}


FUNCTION pois(mean) { 
VERBATIM
    _lpois = - _lmean * log(((double)random()+1.) / ((double)MAXLONG+1.));
ENDVERBATIM
}

FUNCTION poisint(mean) {
  poisint = poisrand(mean)
}

VERBATIM






#define N	16
#define MASK	((unsigned)(1 << (N - 1)) + (1 << (N - 1)) - 1)
#define LOW(x)	((unsigned)(x) & MASK)
#define HIGH(x)	LOW((x) >> N)
#define MUL(x, y, z)	{ long l = (long)(x) * (long)(y); \
		(z)[0] = LOW(l); (z)[1] = HIGH(l); }
#define CARRY(x, y)	((long)(x) + (long)(y) > MASK)
#define ADDEQU(x, y, z)	(z = CARRY(x, (y)), x = LOW(x + (y)))
#define X0	0x330E
#define X1	0xABCD
#define X2	0x1234
#define A0	0xE66D
#define A1	0xDEEC
#define A2	0x5
#define C	0xB
#define SET3(x, x0, x1, x2)	((x)[0] = (x0), (x)[1] = (x1), (x)[2] = (x2))
#define SEED(x0, x1, x2) (SET3(x, x0, x1, x2), SET3(a, A0, A1, A2), c = C)

static unsigned x[3] = { X0, X1, X2 }, a[3] = { A0, A1, A2 }, c = C;
static unsigned short lastx[3];
static void next();

double
my_drand48()
{
	static double two16m = 1.0 / (1L << N);

	next();
	return (two16m * (two16m * (two16m * x[0] + x[1]) + x[2]));
}

static void
next()
{
	unsigned p[2], q[2], r[2], carry0, carry1;

	MUL(a[0], x[0], p);
	ADDEQU(p[0], c, carry0);
	ADDEQU(p[1], carry0, carry1);
	MUL(a[0], x[1], q);
	ADDEQU(p[1], q[0], carry0);
	MUL(a[1], x[0], r);
	x[2] = LOW(carry0 + carry1 + CARRY(p[1], r[0]) + q[1] + r[1] +
		a[0] * x[2] + a[1] * x[1] + a[2] * x[0]);
	x[1] = LOW(p[1] + r[0]);
	x[0] = LOW(p[0]);
}

void
my_srand48(seedval)
long seedval;
{
	SEED(X0, LOW(seedval), HIGH(seedval));
}

#if 0
#ifdef DRIVER

#include <stdio.h>

main()
{
	int i;

	for (i = 0; i < 80; i++) {
		printf("%4d ", (int)(4096 * my_drand48()));
		printf("%.4X%.4X%.4X\n", x[2], x[1], x[0]);
	}
}
#endif
#endif
ENDVERBATIM