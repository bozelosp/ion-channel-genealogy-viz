NEURON {
  SUFFIX stats
  GLOBAL  INSTALLED,seed,kmeasure,verbose,self_ok_combi,hretval,flag,transpose,newline
}

PARAMETER {
  
  INSTALLED=0
  kmeasure=0
  verbose=0
  self_ok_combi=0
  hretval=0
  transpose=0
  newline=90
  flag=0 
}

ASSIGNED { seed }

VERBATIM
#include "misc.h"
#define MIN_MERGESORT_LIST_SIZE    32

union dblint {
  int i[2];
  double d;
};

static u_int32_t ilow=0;  
static u_int32_t ihigh=0; 
static double *x1x, *y1y, *z1z;
static void vprpr();

static int compare_ul(const void* l1, const void* l2) {
  int retval;
  unsigned long d;
  d = (*((const unsigned long*) l1)) - (*((const unsigned long*) l2));
  if(d==0) return 1;
  if(d < 0) return -1;
  return 0;
}

ENDVERBATIM
 


VERBATIM
static double slope(void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2;
	
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) { 
          timestep = *getarg(1); 
        } else { printf("You must supply a timestep\n"); return 0; }

        sigxy= sigx= sigy= sigx2=0; 

        x = (double *) malloc(sizeof(double)*n);
        for(i=0; i<n; i++) {
          x[i] = timestep*i;
          sigxy += x[i] * y[i];
          sigx  += x[i];
          sigy  += y[i];
          sigx2 += x[i]*x[i];
        }
        free(x);
        return (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
}
ENDVERBATIM


VERBATIM
static double moment (void* vv) {
  int i, j, n, fl;
  double *mdata, *y;
  double ave,adev,sdev,svar,skew,curt,s,p;
  n = vector_instance_px(vv, &mdata);
  fl=0;
  if (n<=1) {printf("n must be at least 2 in stats:moment()"); hxe();}
  if(ifarg(1)) {
    if (hoc_is_object_arg(1)) {
      y=vector_newsize(vector_arg(1), 6); fl=1;
    } else { 
      printf("vec.moment(ovec) stores in ovec: ave,adev,sdev,svar,skew,kurt\n");
      return 0;
    }
  }
  for (j=0,s=0;j<n;j++) s+=mdata[j];
  ave=s/n; adev=svar=skew=curt=0.0;
  for (j=0;j<n;j++) { adev+=fabs(s=mdata[j]-ave); svar+=(p=s*s); skew+=(p*=s); curt+=(p*=s); }
  adev/=n;  svar/=(n-1);  sdev=sqrt(svar);
  if (svar) {
    skew /= (n*svar*sdev);
    curt= curt/(n*svar*svar)-3.0;
  } else {printf("No skew/kurtosis when variance = 0 (in stats::moment())\n"); hxe();}
  if (fl) {y[0]=ave; y[1]=adev; y[2]=sdev; y[3]=svar; y[4]=skew; y[5]=curt;}
  return curt;
}
ENDVERBATIM


VERBATIM
static double vslope (void* vv) {
	int i, n;
	double *x, *y;
        double timestep, sigxy, sigx, sigy, sigx2;
	
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) {
          if(vector_arg_px(1, &x) != n ) {
            hoc_execerror("Vector size doesn't match.", 0); 
          }
          sigxy= sigx= sigy= sigx2=0; 

          for(i=0; i<n; i++) {
            sigxy += x[i] * y[i];
            sigx  += x[i];
            sigy  += y[i];
            sigx2 += x[i]*x[i];
          }
        }         
        return (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
}
ENDVERBATIM

VERBATIM







double getsqerr(double* x,double* y,double m,double b,int n,double* meansqerr,double* maxsqerr){
  int i; double val;
  if(!n){
    return -1.0;
  }
  val=0.0;
  *meansqerr=0.0;
  *maxsqerr=0.0;
  for(i=0;i<n;i++){
    val = y[i] - (m*x[i]+b);
    val = val*val;
    if(val>*maxsqerr) *maxsqerr = val;
    *meansqerr += val;
  }
  *meansqerr=*meansqerr/(double)n;
  return *meansqerr;
}
ENDVERBATIM
 

VERBATIM
static double stats(void* vv) {
	int i, n;
	double *x, *y, *out;
        double timestep, sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b, dmeansqerr,dmaxsqerr;
	
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) { 
          timestep = *getarg(1); 
        } else { printf("You must supply a timestep\n"); return 0; }

        sigxy= sigx= sigy= sigx2=sigy2= 0; 

        x = (double *) malloc(sizeof(double)*n);
        for(i=0; i<n; i++) {
          x[i] = timestep*i;
          sigxy += x[i] * y[i];
          sigx  += x[i];
          sigy  += y[i];
          sigx2 += x[i]*x[i];
          sigy2 += y[i]*y[i];
        }
        m = (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
        b = (sigy*sigx2 - sigx*sigxy)/(n*sigx2 - sigx*sigx);
        r = (n*sigxy - sigx*sigy)/(sqrt(n*sigx2-sigx*sigx) * sqrt(n*sigy2-sigy*sigy));
        getsqerr(x,y,m,b,n,&dmeansqerr,&dmaxsqerr); 
        if(ifarg(2)){ 
          out=vector_newsize(vector_arg(2),5);
          out[0]=m; out[1]=b; out[2]=r; out[3]=dmeansqerr; out[4]=dmaxsqerr;
        } else {
          printf("Examined %d data points\n", n);
          printf("slope     = %f\n", m);
          printf("intercept = %f\n", b);
          printf("R         = %f\n", r);
          printf("R-squared = %f\n", r*r);
          printf("MeanSQErr = %f\n",dmeansqerr);
          printf("MaxSQErr  = %f\n",dmaxsqerr);
        }
        free(x);
        return 1;
}

typedef struct pcorst_ {
  int pidse[2];
  double* X;
  double* Y;
  double sigx;
  double sigy;
  double sigx2;
  double sigy2;
  double sigxy;
} pcorst;

void* PCorrelTHFunc(void *arg) {
  pcorst* p;
  int i;
  double *X,*Y;
  p=(pcorst*)arg;

  X = p->X; Y = p->Y;
  p->sigx=p->sigy=p->sigxy=p->sigx2=p->sigy2=0.0;
  for(i=p->pidse[0]; i<p->pidse[1]; i++) {






    p->sigxy += X[i] * Y[i];
    p->sigx += X[i];
    p->sigy += Y[i];
    p->sigx2 += X[i] * X[i];
    p->sigy2 += Y[i] * Y[i];
  }  
  return NULL;
}


#if defined(t)
static double pcorrelsmt(double *x, double* y, int n,int nth) {
  int i,nperth,idx,rc;
  double sigxy, sigx, sigy, sigx2, sigy2, ret;
  pcorst** pp;
  pthread_t* pth;
  pthread_attr_t attr;
  ret=sigxy=sigx=sigy=sigx2=sigy2=0.0; 
  nperth = n / nth;
  
  pp = (pcorst**)malloc(sizeof(pcorst*)*nth);
  idx=0; 
  for(i=0;i<nth;i++) {
    pp[i] = (pcorst*)calloc(1,sizeof(pcorst));
    pp[i]->X = x;
    pp[i]->Y = y;    
    pp[i]->pidse[0] = idx;
    pp[i]->pidse[1] = idx + nperth;
    idx += nperth;
  }
  i--;  if(pp[i]->pidse[1] < n ||
           pp[i]->pidse[1] > n) pp[i]->pidse[1] = n; 
  
  pth=(pthread_t*)malloc(sizeof(pthread_t)*nth);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  for(i=0;i<nth;i++) if((rc=pthread_create(&pth[i], NULL, PCorrelTHFunc, (void*)pp[i]))) {
    printf("pcorrelsmt ERRA: couldn't create thread : %d!\n",rc);
    goto PCMTFREE;
  }
  pthread_attr_destroy(&attr);
  
  for(i=0;i<nth;i++) if((rc=pthread_join(pth[i], NULL))) {
    printf("pcorrelsmt ERRB: couldn't join thread : %d!\n",rc);
    goto PCMTFREE;
  }
  
  for(i=0;i<nth;i++) {
    sigx += pp[i]->sigx;
    sigy += pp[i]->sigy;
    sigxy += pp[i]->sigxy;
    sigx2 += pp[i]->sigx2;
    sigy2 += pp[i]->sigy2;
  }
  sigxy -= (sigx * sigy) / n;
  sigx2 -= (sigx * sigx) / n;
  sigy2 -= (sigy * sigy) / n;
  if(sigx2 <= 0) goto PCMTFREE;
  if(sigy2 <= 0) goto PCMTFREE;
  ret = sigxy / sqrt(sigx2*sigy2);
PCMTFREE:
  
  for(i=0;i<nth;i++) free(pp[i]);
  free(pp);
  free(pth);
  return ret; 
}
#endif


static double pcorrels2 (double *x, double* y, int n) {
  int i;
  double sigxy, sigx, sigy, sigx2, sigy2;
  sigxy=sigx=sigy=sigx2=sigy2=0.0; 
  for(i=0; i<n; i++) {
    sigxy += x[i] * y[i];
    sigx  += x[i];
    sigy  += y[i];
    sigx2 += x[i]*x[i];
    sigy2 += y[i]*y[i];
  }
  sigxy -= (sigx * sigy) / n;
  sigx2 -= (sigx * sigx) / n;
  sigy2 -= (sigy * sigy) / n;
  if(sigx2 <= 0) return 0;
  if(sigy2 <= 0) return 0;
  sigxy = sigxy / sqrt(sigx2*sigy2);
  return sigxy;
}

static double pcorrel (void* vv) {
 int i, n;
  double *x, *y;
  n = vector_instance_px(vv, &x);
  if ((i=vector_arg_px(1, &y)) != n ) {printf("pcorrelsERRA: %d %d\n",n,i); hxe();}
  if(ifarg(2)) {
#if defined(t)
    return pcorrelsmt(x,y,n,(int)*getarg(2));
#else
    printf("using NEURON version 6; pcorrelsmt() not compiled\n"); 
    return 0.0;
#endif
  } else {
    return pcorrels2(x,y,n);
  }
}

ENDVERBATIM







FUNCTION rpval () {
  VERBATIM
  double n , r, df , TINY , ts , mpval;
  n = *getarg(1);
  r = *getarg(2);
  if( r < -1.0 || r > 1.0 ){
    printf("ppval ERRA: r=%g must be : -1.0 <= r <= 1.0\n",r);
    return -1.0;
  }
  if( n < 3 ){
    printf("ppval ERRB: n too small, can't calc probability on samples with < 3 values!\n");
    return -1.0;
  }
  df = n-2; 
  
  
  
  TINY = 1.0e-20;
  ts = r*sqrt(df/((1.0-r+TINY)*(1.0+r+TINY)));
  mpval = betai(0.5*df,0.5,df/(df+ts*ts));
  return mpval;
  ENDVERBATIM
}

VERBATIM


static const double* sortdata = NULL; 


static
int compare(const void* a, const void* b)
{ const int i1 = *(const int*)a;
  const int i2 = *(const int*)b;
  const double term1 = sortdata[i1];
  const double term2 = sortdata[i2];
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}

void csort (int n, const double mdata[], int index[])

{ int i;
  sortdata = mdata;
  for (i = 0; i < n; i++) index[i] = i;
  qsort(index, n, sizeof(int), compare);
}

static double* getrank (int n, double mdata[])

{ int i;
  double* rank;
  int* index;
  rank = malloc(n*sizeof(double));
  if (!rank) return NULL;
  index = malloc(n*sizeof(int));
  if (!index)
  { free(rank);
    return NULL;
  }
  
  csort (n, mdata, index);
  
  for (i = 0; i < n; i++) rank[index[i]] = i;
  
  i = 0;
  while (i < n)
  { int m;
    double value = mdata[index[i]];
    int j = i + 1;
    while (j < n && mdata[index[j]] == value) j++;
    m = j - i; 
    value = rank[index[i]] + (m-1)/2.;
    for (j = i; j < i + m; j++) rank[index[j]] = value;
    i += m;
  }
  free (index);
  return rank;
}


static double spearman(int n, double* data1, double* data2)
{ int i;
  int m = 0;
  double* rank1;
  double* rank2;
  double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double avgrank;
  double* tdata1;
  double* tdata2;
  tdata1 = malloc(n*sizeof(double));
  if(!tdata1) return 0.0; 
  tdata2 = malloc(n*sizeof(double));
  if(!tdata2) 
  { free(tdata1);
    return 0.0;
  }
  for (i = 0; i < n; i++)
  { tdata1[m] = data1[i];
    tdata2[m] = data2[i];
    m++;
  }
  if (m==0) return 0;
  rank1 = getrank(m, tdata1);
  free(tdata1);
  if(!rank1) return 0.0; 
  rank2 = getrank(m, tdata2);
  free(tdata2);
  if(!rank2) 
  { free(rank1);
    return 0.0;
  }
  avgrank = 0.5*(m-1); 
  for (i = 0; i < m; i++)
  { const double value1 = rank1[i];
    const double value2 = rank2[i];
    result += value1 * value2;
    denom1 += value1 * value1;
    denom2 += value2 * value2;
  }
  
  free(rank1);
  free(rank2);
  result /= m;
  denom1 /= m;
  denom2 /= m;
  result -= avgrank * avgrank;
  denom1 -= avgrank * avgrank;
  denom2 -= avgrank * avgrank;
  if (denom1 <= 0) return 0; 
  if (denom2 <= 0) return 0; 
  result = result / sqrt(denom1*denom2);
  return result;
}

static double scorrel(void* vv) {
  int i, n;
  double *x, *y;
  n = vector_instance_px(vv, &x);
  if ((i=vector_arg_px(1, &y)) != n ) {printf("scorrelERRA: %d %d\n",n,i); hxe();}
  return spearman(n,x,y);
}





double Erfcc (double x) {
  double	mt,z,ans;
  z=fabs(x);
  mt=1.0/(1.0+0.5*z);
  ans=mt*exp(-z*z-1.26551223+mt*(1.00002368+mt*(0.37409196+mt*(0.09678418+\
      mt*(-0.18628806+mt*(0.27886807+mt*(-1.13520398+mt*(1.48851587+\
      mt*(-0.82215223+mt*0.17087277)))))))));
  return x >= 0.0 ? ans : 2.0-ans;
}



double Rktau (double* x, double* y, int n){
  int i,j; double c,vx,vy,sx,sy,var,z,tau;
  c = vx = vy = 0.0;
  for(i = 0; i < n; i++) {
    for(j = 0; j < i; j++) {
      sx = (x[i] - x[j]);
      sx = ((sx > 0) ? 1 : ((sx == 0)? 0 : -1));
      sy = (y[i] - y[j]);
      sy = ((sy > 0) ? 1 : ((sy == 0)? 0 : -1));
      vx += sx * sx;
      vy += sy * sy;
      c += sx * sy;
    }
  }  
  if(vx>0 && vy>0) {    
    tau = c / sqrt(vx*vy); 
    return tau;
  }
  return 0.;
}



static double kcorrel (void* vv) {
  int i, n;
  double *x, *y, *prob, *i1d, *i2d, *ps, var, z, tau;
  n = vector_instance_px(vv, &x);
  if ((i=vector_arg_px(1, &y)) != n ) {printf("kcorrel ERRA: %d %d\n",n,i); hxe();}
  if(ifarg(2) && *getarg(2)) {
    i1d=dcrset(n*3); i2d=&i1d[n]; ps=&i2d[n]; tau=kcorfast(x,y,i1d,i2d,n,ps);
  } else {
    tau = Rktau(x,y,n); 
  }
  if(!(ifarg(3) && vector_arg_px(3,&prob))) prob = 0x0; 
  if(prob) { 
    var = (4.0 * n + 10.0) / (9.0 * n * (n - 1.0));
    z = tau / sqrt(var);
    *prob = Erfcc(fabs(z)/1.4142136); 
  }
  return tau;
}


static int mycompare (const void* a, const void* b)
{ double d1,d2;
  d1 = *(double*)a;
  d2 = *(double*)b;
  if (d1 < d2) return -1;
  if (d1 > d2) return +1;
  return 0;
}



void mergesort_array (double a[], int size, double temp[],unsigned long* swapcount) {
  int i1, i2, i, tempi, j, vv;
  double *right,*left;
  if(size<=1) return; 
  if (0 && size <MIN_MERGESORT_LIST_SIZE){
    
    for (i=0; i < size; i++) {
      vv = a[i];
      for (j = i - 1; j >= 0; j--) {
        if (a[j] <= vv) break;
        a[j + 1] = a[j];
      }
      a[j + 1] = vv;
    }
    return;
  }  
  mergesort_array(a, size/2, temp,swapcount);                 
  mergesort_array(a + size/2, size - size/2, temp,swapcount); 
  
  i=tempi=0;  i1 = size/2; i2 = size - size/2;
  left = a; right = &a[size/2];
  while(i1>0 && i2>0) { 
    if(*right < *left) {
      *swapcount += i1; 
      temp[i] = *right++;
      i2--;
    } else {
      temp[i] = *left++;
      i1--;
    }
    i++;
  }
  if(i2>0) { 
    while(i2-->=0 && i<size) temp[i++] = *right++; 
  } else {
    while(i1-->=0 && i<size) temp[i++] = *left++; 
  }
  memcpy(a, temp, size*sizeof(double));
}




int qsort2 (double *p1in, double* p2in, int n,double* p1out,double* p2out) {
  int i;
  scr=scrset(n);
  for (i=0;i<n;i++) scr[i]=i;
  nrn_mlh_gsort(p1in, scr, n, cmpdfn);
  for (i=0;i<n;i++) {
    p1out[i]=p1in[scr[i]];
    p2out[i]=p2in[scr[i]];
  }
  return 1;
}


unsigned long getMs (double* data,int n) {  
  unsigned long Ms, tieCount;
  int i;
  Ms = tieCount = 0;
  for(i=1;i<n;i++) {
    if(data[i] == data[i-1]) {
      tieCount++;
    } else if(tieCount) {
      Ms += (tieCount*(tieCount+1))/2;
      tieCount = 0;
    }
  }
  if(tieCount) {
    Ms += (tieCount*(tieCount+1)) / 2;
  }
  return Ms;
}





double kcorfast (double* input1, double* input2, double* i1d , double* i2d,int n,double* ps) {
    int i;
    unsigned long nPair, N, m1, m2, tieCount, swapCount;
    long s;
    double denom1,denom2;
    m1 = m2 = 0; N = n;
    nPair = N * ( N - 1 ) / 2; 
    qsort2(input1,input2,n,i1d,i2d); 
    s = nPair; 

    if(verbose>2) printf("nPair=%lu\n",nPair);
    if(verbose>3){printf("i1d after qsort2: "); for(i=0;i<n;i++) printf("%g ",i1d[i]); printf("\n");
                  printf("i2d after qsort2: "); for(i=0;i<n;i++) printf("%g ",i2d[i]); printf("\n");}
    tieCount = 0;
    for(i=1;i<n;i++) {
        if(i1d[i] == i1d[i-1]) {
            tieCount++;
        } else if(tieCount > 0) {
            qsort(&i2d[i-tieCount-1],tieCount+1,sizeof(double),mycompare);
            m1 += tieCount * (tieCount + 1) / 2;
            s += getMs(&i2d[i-tieCount-1],tieCount+1);
            tieCount = 0;
        }
    }
    if(verbose>2) printf("tieCount=%lu\n",tieCount);
    if(tieCount > 0) {
        qsort(&i2d[n-tieCount-1],tieCount+1,sizeof(double),mycompare);
        m1 += tieCount * (tieCount + 1) / 2;
        s += getMs(&i2d[n-tieCount-1],tieCount+1);
    }
    if(verbose>2) printf("tieCount=%lu\n",tieCount);
    swapCount = 0;

    mergesort_array(i2d,n,ps,&swapCount); 
    if(verbose>3) { printf("i2d after mergesort: "); for(i=0;i<n;i++) printf("%g ",ps[i]); printf("\n"); }
    if(verbose>2) printf("swapCount=%lu\n",swapCount);

    m2 = getMs(i2d,n); if(verbose>2) printf("s=%lu m1=%lu m2=%lu\n",s,m1,m2);
    s -= (m1 + m2) + 2 * swapCount; 
    denom1=nPair-m1; denom2=nPair-m2; if(verbose>2) printf("s=%lu d1=%g d2=%g\n",s,denom1,denom2);
    if(denom1>0. && denom2>0.) return s / sqrt(denom1*denom2); else return 0.;
}

static double rms (void* vv) {
  int i,n;
  double *x,sum;
  if(!(n=vector_instance_px(vv, &x))) {printf("rms ERRA: 0 sized vector!\n"); hxe();}
  sum=0.0;
  for(i=0;i<n;i++) sum += x[i]*x[i];
  sum/=(double)n;
  if(sum>0.) return sqrt(sum); else return 0.0;
}


static double cumsum (void* vv) {
  int i,n;
  double *x,*y;
  if(!(n=vector_instance_px(vv, &x))) {printf("cumsum ERRA: 0 sized vector!\n"); hxe();}
  if(vector_arg_px(1, &y) != n) {printf("cumsum ERRB: output vec size needs size of %d\n",n); hxe();}
  memcpy(y,x,sizeof(double)*n);
  for(i=1;i<n;i++) y[i] += y[i-1];
  return 1.0;
}

ENDVERBATIM
 

VERBATIM
static double unnan (void *vv) {
  int i,nx,cnt; double newnan,newinf,neginf;
  union dblint xx;
  double *x;
  newnan=newinf=neginf=0;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) newinf=newnan=*getarg(1);
  if (ifarg(2)) newinf=*getarg(2);
  if (ifarg(3)) neginf=*getarg(3);
  for (i=0,cnt=0;i<nx;i++) { 
    xx.d=x[i];
    if (xx.i[0]==0x0 && xx.i[1]==0xfff80000) {x[i]=newnan; cnt++;}
    if (xx.i[0]==0x0 && xx.i[1]==0x7ff00000) {x[i]=newinf; cnt++;}
    if (xx.i[0]==0x0 && xx.i[1]==0xfff00000) {x[i]=neginf; cnt++;}
  }
  return (double)cnt;
}
ENDVERBATIM


VERBATIM
static double vstats(void* vv) {
	int i, n;
	double *x, *y, *out;
        double sigxy, sigx, sigy, sigx2, sigy2;
        double r, m, b, dmeansqerr,dmaxsqerr;
	
	n = vector_instance_px(vv, &y);

        if(ifarg(1)) {
          if(vector_arg_px(1, &x) != n ) {
            hoc_execerror("Vector size doesn't match.", 0); 
          }
          sigxy= sigx= sigy= sigx2=sigy2=0; 

          for(i=0; i<n; i++) {
            sigxy += x[i] * y[i];
            sigx  += x[i];
            sigy  += y[i];
            sigx2 += x[i]*x[i];
            sigy2 += y[i]*y[i];
          }
          m = (n*sigxy - sigx*sigy)/(n*sigx2 - sigx*sigx);
          b = (sigy*sigx2 - sigx*sigxy)/(n*sigx2 - sigx*sigx);
          r = (n*sigxy - sigx*sigy)/(sqrt(n*sigx2-sigx*sigx) * sqrt(n*sigy2-sigy*sigy));
          getsqerr(x,y,m,b,n,&dmeansqerr,&dmaxsqerr);
          if(ifarg(2)){ 
            out=vector_newsize(vector_arg(2),5);
            out[0]=m; out[1]=b; out[2]=r; out[3]=dmeansqerr; out[4]=dmaxsqerr;
          } else {
            printf("Examined %d data points\n", n);
            printf("slope     = %f\n", m);
            printf("intercept = %f\n", b);
            printf("R         = %f\n", r);
            printf("R-squared = %f\n", r*r);
            printf("MeanSQErr = %f\n",dmeansqerr);
            printf("MaxSQErr  = %f\n",dmaxsqerr);
          }
          return 1;
        } else {
          printf("You must supply an x vector\n");
          return 0;
        }
}
ENDVERBATIM



VERBATIM
static double randwd(void* vv) {
	int i, ii, jj, nx, ny, flip, flag;
	double* x, *y;
	
	nx = vector_instance_px(vv, &x);
        flip = (int) *getarg(1);
        if (ifarg(2)) { 
          flag = 1; ny = vector_arg_px(2, &y);
          if (ny!=flip) { hoc_execerror("Opt vector must be size for # of flips", 0); }
        } else { flag = 0; }
        if (flip>=nx) { hoc_execerror("# of flips exceeds (or ==) vector size", 0); }
	for (i=0; i < nx; i++) { x[i] = BVBASE; }
	for (i=0,jj=0; i < flip; i++) { 
	  ii = (int) ((nx+1)*drand48());
	  if (x[ii]==BVBASE) {
	    x[ii] = 1.; 
            if (flag) { y[jj] = ii; jj++; }
	  } else {
	    i--;
	  }
	}
	return flip;
}
ENDVERBATIM
 

VERBATIM
static double hash (void* vv) {
  int i, j, nx, nv[VRRY], num, vfl; 
  union dblint xx;
  Object* ob;
  double *x, *vvo[VRRY], big, y, prod;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) {
    vfl=0;
    ob=*hoc_objgetarg(1); 
    num = ivoc_list_count(ob);
    if (num>VRRY) {printf("vecst:hash ERR: can only handle %d vecs: %d\n",VRRY,num); hxe();}
    for (i=0;i<num;i++) { nv[i] = list_vector_px(ob, i, &vvo[i]);
      if (nx!=nv[i]) { printf("vecst:hash ERR %d %d %d\n",i,nx,nv[i]);hxe();}}
  } else {
    vfl=1; num=nx; nx=1; 
  }
  big=pow(DBL_MAX,1./(double)num); 
  for (i=0;i<nx;i++) {
    for (j=0,prod=1;j<num;j++) {
      if (vfl) {  xx.d=x[j];       
      } else   {  xx.d=vvo[j][i]; }
      if (xx.i[0]==0) { xx.i[0]=xx.i[1]; xx.i[0]<<=4; } 
      if (xx.i[1]==0) { xx.i[1]=xx.i[0]; xx.i[1]<<=4; } 
      mcell_ran4_init(&xx.i[1]);
      mcell_ran4(&xx.i[0], &y, 1, big); 
      prod*=y;  
    }
    if (! vfl) x[i]=prod; else return prod; 
  }
  return (double)nx;
}
ENDVERBATIM





VERBATIM
static double smash (void* vv) {
  int i, j, nx, nv[VRRY], num; 
  Object* ob;
  double *x, *vvo[VRRY], wt, wtj;
  nx = vector_instance_px(vv, &x);
  ob=*hoc_objgetarg(1); 
  if (ifarg(2)) wtj=*getarg(2); else wtj=10.;
  num = ivoc_list_count(ob);
  if (num>VRRY) {printf("vecst:smash ERRA: can only handle %d vecs: %d\n",VRRY,num); hxe();}
  if (transpose) if (nx!=num) { printf("vecst:smash ERRB %d %d %d\n",i,nx,nv[i]);hxe(); }
  for (i=0;i<num;i++) { 
    nv[i] = list_vector_px(ob, i, &vvo[i]);
    if (!transpose) if (nx!=nv[i]) { printf("vecst:smash ERRB %d %d %d\n",i,nx,nv[i]);hxe(); }
  }
  if (transpose) { 
    for (i=0;i<num;i++) { 
      for (j=0,x[i]=0,wt=1;j<nv[i];j++,wt*=wtj) x[i]+=vvo[i][j]*wt;
    }
  } else for (i=0;i<nx;i++) {
    for (j=0,x[i]=0,wt=1;j<num;j++,wt*=wtj) x[i]+=vvo[j][i]*wt;
  }
  return (double)nx;
}
ENDVERBATIM





VERBATIM
static double smash1 (void* vv) {
  int i, j, nx, nv[VRRY], num, mod; 
  double *x, wt, wtj, res;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) wtj=*getarg(1); else wtj=10.;
  if (ifarg(2)) mod=(int)*getarg(2); else mod=0;
  for (j=0,res=0,wt=1;j<nx;j++,wt*=wtj) {
    res+=x[j]*wt;
    if (mod && j%mod==0) wt=1;
  }
  return res;
}
ENDVERBATIM


VERBATIM 
static double dpro (void* vv) {
  int i, j, nx, nv[VRRY], num, step, gap; 
  Object* ob;
  double *x, *vvo[VRRY], wt;
  nx = vector_instance_px(vv, &x);
  ob=*hoc_objgetarg(1); 
  if (ifarg(2)) step=(int)*getarg(2); else step=1;
  if (ifarg(3)) gap=(int)*getarg(3); else gap=1;
  num = ivoc_list_count(ob);
  if (num>VRRY) {printf("stats:dpro ERR: can only handle %d vecs: %d\n",VRRY,num); hxe();}
  for (i=0;i<num;i++) { 
    nv[i] = list_vector_px(ob, i, &vvo[i]);
    if (nx!=nv[i]) { printf("stats:dpro ERR %d %d %d\n",i,nx,nv[i]);hxe(); }
  }
  for (i=0;i<nx;i+=step) {
    for (j=0,x[i]=0,wt=1;j<num;j++) {
      x[i]+=vvo[j][i]*wt;
    }
  }
  return (double)nx;
}
ENDVERBATIM














VERBATIM
static double setrnd (void* vv) {
  int flag, i,j,k,n,cnt; unsigned int nx, nx1, nex, lt, rt, mid;
  double *x, y, *ex, *ex2, min, max, dfl, tmp, step, num;
  unsigned long value;
  value=1;
  nx = vector_instance_px(vv, &x);
  flag = (int)(dfl=*getarg(1));
  if (flag==1) {
    for (i=0; i < nx; i++) x[i] = (double)rand()/RAND_MAX; 
  }  else if (flag==2) {
    for (i=0; i < nx; i++) x[i] = drand48(); 
  } else if (flag==3) { 
    unsigned long a = 2147437301, c = 453816981, m = ~0;
    value = (unsigned long) seed;
    for (i=0; i < nx; i++) {
      value = a * value + c;
      x[i] = (fabs((double) value / (double) m));
    }
    seed=(double)value;
  } else if (flag==4) { 
    ex=0x0; i=2;
    if (ifarg(i)) {
      if (hoc_is_object_arg(i)) {
        nex=vector_arg_px(i++,&ex); 
        step=ifarg(i)?*getarg(i):1.0;
        max=1.0; i++;
      } else {
	if (dfl==4.5 || ifarg(4)) { 
	  min=*getarg(i++); 
	  max=*getarg(i++)-min; 
	  dfl=4.5;
	} else {
	  max=*getarg(i++);
	}
      }
    } else max=1.0; 
    if (ifarg(i)) { y=*getarg(i++); if (y) ihigh=(unsigned int)y; } 
    if (max==0) { for (i=0;i<nx;i++) x[i]=0.;
    } else mcell_ran4(&ihigh, x, nx, max);
    if (dfl==4.5) for (i=0;i<nx;i++) x[i]+=min;
    if (ex) for (i=0;i<nx;i++) { 
      num=x[i]; lt=0; rt=nex-1;
      while (lt<=rt) { 
        mid=(lt+rt)/2;
        if (num>ex[mid]) { 
          if (num<ex[mid+1]) break; 
          lt=mid+1;
        } else if (num<ex[mid]) rt=mid-1;
      }
      x[i]=step*mid;
    }
    return (double)ihigh;
  } else if (flag==5) { 
    n=100; ex=0x0;
    if (ifarg(2)) {
      if (hoc_is_object_arg(2)) {
        n=vector_arg_px(2,&ex); 
      } else {
        n=(int)*getarg(2);
      }
    }
    i=3; 
    if (dfl==5.5 || ifarg(4)) { max=*getarg(3); min=n; n=max-min+1; dfl=5.5; i=4; }
    if (ifarg(i)) { y=*getarg(i); if (y) ihigh=(unsigned int)y; }
    if (n<=1) { for (i=0;i<nx;i++) x[i]=0.0;
    } else mcell_ran4(&ihigh, x, nx, (double)n);
    if (dfl==5.5)  { for (i=0;i<nx;i++) x[i]=min+floor(x[i]);
    } else if (ex) { for (i=0;i<nx;i++) x[i]=  ex[(int)x[i]];
    } else           for (i=0;i<nx;i++) x[i]=    floor(x[i]);
    return (double)ihigh;
  } else if (flag==6) { 
    min=*getarg(2); max=*getarg(3); i=4; nex=0;
    if (ifarg(i+1)) {
      nex=vector_arg_px(i, &ex); 
      ihigh=(unsigned int)(*getarg(i+1));
    } else if (ifarg(i)) {
      if (hoc_is_object_arg(i)) { nex=vector_arg_px(i, &ex); 
      } else { ihigh=(unsigned int)(*getarg(i)); }
    } 
    
    if (nex>1) { 
      scrset(nex);
      x1x = (double *)realloc(x1x,sizeof(double)*nx*4);
      for (i=0;i<nex;i++) scr[i]=i;
      nrn_mlh_gsort(ex, scr, nex, cmpdfn);
      for (i=0;i<nex;i++) x1x[i]=ex[scr[i]];
      for (i=0;i<nex;i++) ex[i]=x1x[i];
    }
    for (j=0;j<nex;j++) if (ex[j]>max || ex[j]<min || (j>0 && ex[j]<=ex[j-1])) {
      printf("%g in exclusion list -- out of range [%g,%g] or a repeat\n",ex[j],min,max);hxe();} 
    if (max-min+1-nex==nx) { 
      for (i=min,k=0;i<=max;i++) {
        y=(double)i;
        for (j=0;j<nex && ex[j]!=y;j++) {} 
        if (j==nex) x[k++]=y; 
      }
      dshuffle(x,nx);
      return (double)nx;
    } else if (max-min+1-nex<nx) {
      printf("setrnd ERR6A incompatible: min=%g max=%g; %d vs %g\n",min,max,nx,max-min+1-nex); 
      hxe();
    }
    cnt=0; nx1=nx; 
    while (cnt<nx && nx1<=256*nx) {
      nx1*=4; 
      x1x = (double *)realloc(x1x,sizeof(double)*nx1);
      y1y = (double *)realloc(y1y,sizeof(double)*nx1);
      z1z = (double *)realloc(z1z,sizeof(double)*nx1);
      mcell_ran4(&ihigh, x1x, nx1, max-min+1-nex);
      for (i=0;i<nx1;i++) x1x[i]=floor(x1x[i])+min;
      cnt=uniq2(nx1,x1x,y1y,z1z);
    }
    if (nex) { 
      
      for (i=0;i<nx;i++) {
        for (j=0,k=0;j<nex;j++) if (z1z[i]+k>=ex[j]) k++; 
        x[i]=z1z[i]+k;
      }
    } else  for (i=0;i<nx;i++) x[i]=z1z[i];
  }
  return nx;
}
ENDVERBATIM
 



VERBATIM
static double hamming (void* vv) {
  int i, nx, ny, nz, prflag;
  double* x, *y, *z,sum;
  sum = 0.;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  if (ifarg(2)) { 
    prflag = 1; nz = vector_arg_px(2, &z);
  } else { prflag = 0; }
  if (nx!=ny || (prflag && nx!=nz)) {
    hoc_execerror("Vectors must be same size", 0);
  }
  for (i=0; i < nx; ++i) {
    if (x[i] != y[i]) { sum++; 
      if (prflag) { z[i] = 1.; }
    } else if (prflag) { z[i] = 0.; }
  }
  return sum;
}
ENDVERBATIM





VERBATIM
static double combi (void* vv) {
  int i,j,k,m,n,prix,prixe,poix,poixe,prn,pon,s,vfl,tot,soc; 
  double perc, *v1, *v2, *vpr, *vpo, *vec; void *v1v, *v2v;
  int nx,cnt,nx1,nvec,nv1,nv2;
  nx=nvec = vector_instance_px(vv, &vec);
  if (ifarg(7)) vfl=0; else vfl=1; 
  nv1 = vector_arg_px((vfl?3:5), &v1); v1v=vector_arg((vfl?3:5)); 
  nv2 = vector_arg_px((vfl?4:6), &v2); v2v=vector_arg((vfl?4:6)); 
  perc = *getarg(vfl?5:7);
  if (nv1!=nv2) {printf("stats:combi()ERRA out vec size discrep: %d,%d\n",nv1,nv2); hxe();}
  if (vfl) {
    prn = vector_arg_px(1, &vpr);
    pon = vector_arg_px(2, &vpo);
  } else {
    prix=(int)*getarg(1); prixe=(int)*getarg(2); poix=(int)*getarg(3); poixe=(int)*getarg(4);
    prn=prixe-prix+1; pon=poixe-poix+1;
  } 
  if (prn<=0 || pon<=0) {printf("stats:combi()ERRB %d,%d\n",prn,pon); hxe();}
  soc=(int)self_ok_combi;
  if (soc) tot=prn*pon; else { 
    if (vfl){for(i=0,tot=0;   i<prn;   i++) for (j=0;j<pon;j++) if (vpr[i]!=vpo[j]) tot++;
    } else  for (i=prix,tot=0;i<=prixe;i++) for (j=poix;j<=poixe;j++) if (i!=j)     tot++;
  }
  
  if (perc<1) s=(int)floor(perc*(double)tot+0.5); else s=(int)perc; 
  
  if (soc && s==1) { 
    vec=vector_newsize(vv,1); v1=vector_newsize(v1v,nv1+1); v2=vector_newsize(v2v,nv2+1);
    mcell_ran4(&ihigh, vec, 1, prn);
    if (vfl) v1[nv1]=vpr[(int)vec[0]]; else v1[nv1]=prix+floor(vec[0]);
    mcell_ran4(&ihigh, vec, 1, pon);    
    if (vfl) v2[nv2]=vpo[(int)vec[0]]; else v2[nv2]=poix+floor(vec[0]);    
    return 1.0; 
  }
  vec=vector_newsize(vv,s); 
  if (tot==s) { for (i=0;i<s;i++) vec[i]=(double)i; 
  } else { 
    cnt=0; nx1=10*s;
    while (cnt<s && nx1<=640*nx) {
      nx1*=4; 
      x1x = (double *)realloc(x1x,sizeof(double)*nx1);
      y1y = (double *)realloc(y1y,sizeof(double)*nx1);
      z1z = (double *)realloc(z1z,sizeof(double)*nx1);
      mcell_ran4(&ihigh, x1x, nx1, tot);
      for (i=0;i<nx1;i++) x1x[i]=floor(x1x[i]);
      cnt=uniq2(nx1,x1x,y1y,z1z);
    }
    for (i=0;i<s;i++) vec[i]=z1z[i];
  }
  v1=vector_newsize(v1v,nv1+s); v2=vector_newsize(v2v,nv2+s);
  
  for (i=0,m=nv1,n=-1;i<prn;i++) for (j=0;j<pon;j++) { 
    if (vfl) {if (soc || (vpr[i]!=vpo[j])) n++; else continue; 
    } else   {if (soc || (prix+i!=poix+j)) n++; else continue; }
    for (k=0;k<s;k++) if (vec[k]==(double)n) { 
      if (vfl) {v1[m]=vpr[i]; v2[m]=vpo[j];    
      } else   {v1[m]=prix+i; v2[m]=poix+j;   }
      m++;     
      break;
    }
  }
  return (double)s;
}
ENDVERBATIM

VERBATIM


void dshuffle (double* x,int nx) {
  int n,k; double temp,y[1];
  for (n=nx;n>1;) {
    mcell_ran4(&ihigh, y, 1, n);
    n--;
    k=(int)y[0]; 
    temp = x[n];
    x[n] = x[k];
    x[k] = temp;
  }  
}


void uishuffle(unsigned int* x,int nx) {
  int n,k; unsigned int temp; double y[1];
  for (n=nx;n>1;) {
    mcell_ran4(&ihigh, y, 1, n);
    n--;
    k=(int)y[0]; 
    temp = x[n];
    x[n] = x[k];
    x[k] = temp;
  }  
}


void ishuffle(int* x,int nx) {
  int n,k,temp; double y[1];
  for (n=nx;n>1;) {
    mcell_ran4(&ihigh, y, 1, n);
    n--;
    k=(int)y[0]; 
    temp = x[n];
    x[n] = x[k];
    x[k] = temp;
  }  
}

unsigned long choose (int n, int k) {
  int i,delta;
  unsigned long ret;
  
  if (n < k) return 0;
  if (n==k)  return 1;
  if (k < n - k) {
    delta = n - k;
  } else {
    delta = k;
    k = n - k;
  }
  ret = delta + 1;
  for (i = 2; i <= k; ++i)
    ret = (ret * (delta + i)) / i;
  return ret;
}


unsigned long syncci (int nn,int kk,int* ccvv) {
  unsigned long c = 0;
  while ((kk > 0) && (*ccvv >= kk)) {
    c += choose (*ccvv++, kk--);
  }
  return c;
}



unsigned long syncc (int nn,int kk,double* ccvv) {
  unsigned long c = 0;
  while ((kk > 0) && (*ccvv >= kk)) {
    c += choose (*ccvv++, kk--);
  }
  return c;
}


void synccv (int nn,int kk,int cc,double* ccvv) {
  unsigned long n_k;
  while (--nn >= 0) {
    n_k = choose (nn, kk);
    if (cc >= n_k) {
      cc -= n_k;
      *ccvv++ = nn;
      --kk;
    }
  }
}

ENDVERBATIM



FUNCTION combs () {
  VERBATIM
  unsigned int n,k;
  n=(unsigned int)*getarg(1);
  k=(unsigned int)*getarg(2);
  return choose(n,k);
  ENDVERBATIM
}




VERBATIM
static double comb (void* vv) {
  double* x;
  int nn,kk,cc,sz,i;
  sz=vector_instance_px(vv,&x);
  nn=(int)*getarg(1);
  kk=(int)*getarg(2);
  cc=(int)*getarg(3);
  if(sz<kk){
    printf("comb ERRA: output vec sz must be >= %d , is %d!\n",kk,sz);
    return 0.0;
  }
  memset(x,0,sizeof(double)*kk);
  synccv(nn,kk,cc,x);
  vector_resize(vv,kk);
  return 1.0;
}
ENDVERBATIM





VERBATIM
static double combid (void* vv) {
  double* x;
  int nn,kk,sz,i;
  sz=vector_instance_px(vv,&x);
  nn=(int)*getarg(1);
  kk=(int)*getarg(2);
  if(sz<kk){
    printf("comb ERRA: input vec sz must be >= %d , is %d!\n",kk,sz);
    return 0.0;
  }
  return syncc(nn,kk,x);
}
ENDVERBATIM

VERBATIM
int findlong (unsigned long* p,unsigned long val,int istart,int iend) {
  int i;
  for(i=istart;i<=iend;i++) if(p[i]==val) return 1; 
  return 0;
}

ENDVERBATIM













VERBATIM
static double rsampsig(void* vv){
  int n0,n1,na,nn,kk,cc,i,j,*pm,szthis,onesided,nocmbchk,bti,*pids;
  unsigned long nruncombs,nallcombs,*pcombids;
  double *g0,*g1,*ga,prc,*g0t,*g1t,dmobs,dm0,dm1,*phso,nmatch,*pthis,dret;
  void* vhso; 
  Symbol* pHocVecFunc,*pHocCompFunc; 
  dret=-1.0;
  g0t=g1t=NULL; pm=pids=NULL; pcombids=NULL;
  szthis=vector_instance_px(vv, &pthis); 
  n0=vector_arg_px(1,&g0);
  n1=vector_arg_px(2,&g1);
  na=n0+n1;
  prc=*getarg(3);
  if(!(pHocVecFunc=hoc_lookup(gargstr(4)))){
    printf("rsampsig ERRA: couldn't find hoc vec func %s\n",gargstr(4));
    goto CLEANRSAMPSIG;
  }
  if(!(pHocCompFunc=hoc_lookup(gargstr(5)))){
    printf("rsampsig ERRB: couldn't find hoc comp func %s\n",gargstr(5));
    goto CLEANRSAMPSIG;
  }
  if(vector_arg_px(6,&phso)<na){ 
    printf("rsampsig ERRC: arg 6 must have size >= %d!\n",na);
    goto CLEANRSAMPSIG;
  }
  if(prc<=0.0) {
    printf("rsampsig ERRD: invalid value for arg 3, must be > 0.0!\n");
    goto CLEANRSAMPSIG;
  }
  vhso=vector_arg(6);
  ga=(double*)malloc(sizeof(double)*na);
  memcpy(ga,g0,sizeof(double)*n0);
  memcpy(ga+n0,g1,sizeof(double)*n1);
  
  g0t=(double*)malloc(sizeof(double)*n0);
  g1t=(double*)malloc(sizeof(double)*n1);
  if(n0<n1) kk=n0; else kk=n1;
  if(verbose>1) printf("choose(%d,%d)=%ld\n",na,kk,choose(na,kk));
  nallcombs=choose(na,kk);
  nruncombs=prc>1.0?prc:prc*nallcombs; nmatch=0.0;
  if(szthis<nruncombs){
    printf("rsampsig ERRE: vector size (%d) < nruncombs (%ld)!\n",szthis,nruncombs);
    goto CLEANRSAMPSIG;
  }
  onesided=ifarg(7)?(int)*getarg(7):1;
  nocmbchk=ifarg(8)?(int)*getarg(8):1;
  pids=(int*)malloc(sizeof(int)*na); 
  for(i=0;i<na;i++) pids[i]=i; 
  pcombids=(unsigned long*)malloc(sizeof(unsigned long)*nruncombs); 
  if(verbose>1) printf("na=%d , kk=%d, n0=%d, n1=%d\n",na,kk,n0,n1);
  if(verbose>1) printf("nruncombs = %ld\n",nruncombs);
  if(verbose>2) pm=(int*)malloc(sizeof(int)*na);
  for(i=0;i<nruncombs;i++) {    
    do { 
      ishuffle(pids,na);
      pcombids[i] = syncci(na , kk, pids);
    } while(!nocmbchk && i-1>0 && findlong(pcombids,pcombids[i],0,i-1)); 
    if(verbose) if(i%100==0) { printf("."); fflush(stdout); }
    for(j=0;j<n0;j++) *g0t++=ga[*pids++]; 
    for(;j<na;j++)    *g1t++=ga[*pids++]; 
    g0t-=n0; g1t-=n1; pids-=na; 
    if(verbose>2){ for(j=0;j<n0;j++) pm[pids[j]]=0; for(;j<na;j++) pm[pids[j]]=1;  printf("pm: ");
      for(j=0;j<40;j++) printf("%d",pm[j]); printf("\n"); }
    
    vector_resize(vhso, n0); memcpy(phso,g0t,sizeof(double)*n0); 
    hoc_call_func(pHocVecFunc,0); dm0 = hretval; 
    vector_resize(vhso, n1); memcpy(phso,g1t,sizeof(double)*n1); 
    hoc_call_func(pHocVecFunc,0); dm1 = hretval; 
    hoc_pushx(dm0); hoc_pushx(dm1); hoc_call_func(pHocCompFunc,2); 
    pthis[i]=onesided?hretval:fabs(hretval); 
  }
  vector_resize(vv,nruncombs);
  
  vector_resize(vhso,n0); memcpy(phso,g0,sizeof(double)*n0); 
  hoc_call_func(pHocVecFunc,0); dm0 = hretval; 
  vector_resize(vhso,n1); memcpy(phso,g1,sizeof(double)*n1); 
  hoc_call_func(pHocVecFunc,0); dm1 = hretval; 
  hoc_pushx(dm0); hoc_pushx(dm1); hoc_call_func(pHocCompFunc,2); 
  dmobs = onesided?hretval:fabs(hretval); 
  
  for(i=0;i<nruncombs;i++) if(pthis[i] > dmobs) nmatch++;    
  dret=nmatch/(double)nruncombs;
CLEANRSAMPSIG: 
  if(ga) free(ga); if(g0t) free(g0t); if(g1t) free(g1t); if(pm) free(pm);
  if(pcombids) free(pcombids); if(pids) free(pids);
  return dret;
}
ENDVERBATIM


VERBATIM






static double rantran (void* vv) {
  int i,j,ix,ixe,ixvn,nvn,rvn,na,xj;
  double *ixv, *nv, *x, y[1], ixn,step,indx;
  rvn=vector_instance_px(vv, &x);
  for (na=1;ifarg(na);na++); na--; 
  for (i=1;i<na;i+=2) {
    if (hoc_is_object_arg(i)) {
      step=-1;
      ixvn=vector_arg_px(i,&ixv);
      nvn=vector_arg_px(i+1, &nv); 
    } else { 
      step=*getarg(i);
      indx=*getarg(i+1);
      if (indx>1.) { 
        x1x = (double *)realloc(x1x,sizeof(double)*rvn); 
        mcell_ran4(&ihigh, x1x, rvn, indx);
      }
    }
    for (j=0;j<rvn;j++) { 
      if (step>-1) { 
        x[j]=step + x[j]*indx + ((indx>1.)?(floor(x1x[j])):0.0);
      } else {
        xj=(int)x[j]; 
        ix=(int)ixv[xj]; ixn=ixv[xj+1]-ix; 
        if (ixn==1.) {
          x[j]=nv[ix];
        } else {
          mcell_ran4(&ihigh, y, 1, ixn);  
          x[j]=nv[(int)y[0]+ix];
        }
      }
    }      
  }
  return (double)rvn; 
}
ENDVERBATIM




VERBATIM
static double shuffle (void* vv) {
  int i,j,k,n,nx,augfac; double *x, y[1], temp, augstep;
  nx=vector_instance_px(vv, &x);
  if (ifarg(1)) {
    augfac=(int)*getarg(1);
    if (ifarg(2)) augstep=*getarg(2); else augstep=1.0/augfac;
    x=vector_newsize(vv,nx*augfac);
    for (i=1;i<augfac;i++) for (j=0;j<nx;j++) x[i*nx+j]=x[j]+i*augstep;
    nx*=augfac;
  }
  dshuffle(x,nx);
  return (double)nx;
}
ENDVERBATIM


VERBATIM
static double distance (void* vv) {
  int i, nx, ny;
  double* x, *y, sum;
  sum = 0.;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  if (nx!=ny) {printf("Vectors must be same size %d %d\n",nx,ny); hxe();}
  for (i=0; i<nx; i++) sum+=(x[i]-y[i])*(x[i]-y[i]); 
  return sqrt(sum);
}
ENDVERBATIM


VERBATIM
static double ndprd (void* vv) {
  int i, nx, ny;
  double* x, *y, sum, sumx, sumy;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  if (nx!=ny) {printf("Vectors must be same size %d %d\n",nx,ny); hxe();}
  for (i=0, sum=0., sumx=0., sumy=0.; i<nx; i++) {
    sum+=x[i]*y[i]; sumx+=x[i]*x[i]; sumy+=y[i]*y[i]; 
  }
  if (ifarg(2)) { return sum/sqrt(sumx)/sqrt(sumy);                   
  } else {        return acos(sum/sqrt(sumx)/sqrt(sumy))*180./M_PI; } 
}
ENDVERBATIM



VERBATIM
static double flipbits (void* vv) {	
  int i, j, nx, ny, flip, ii;
  double *x, *y;
  nx = vector_instance_px(vv, &x);
  ny = vector_arg_px(1, &y);
  flip = (int)*getarg(2);
  if (ny<nx) {hoc_execerror("flipbits:Scratch vector must adequate size", 0);}
  mcell_ran4(&ihigh, y, (unsigned int)ny, (double)nx); 
  for (i=0,j=0; i<flip && j<ny; j++) { 
    ii=(int)y[j];
    if        (x[ii]==BVBASE) { x[ii]= 1e9; i++;  
    } else if (x[ii]==1)      { x[ii]=-1e9; i++; }
  }
  j=i;
  for (i=0; i<nx; i++) if (x[i]==1e9) x[i]=1; else if (x[i]==-1e9) x[i]=BVBASE;
  return (double)j;
}
ENDVERBATIM
 



VERBATIM
static double flipbalbits(void* vv) {	
	int i, nx, ny, flip, ii, next;
	double* x, *y;

	nx = vector_instance_px(vv, &x);
	ny = vector_arg_px(1, &y);
        flip = (int)*getarg(2);
	if (nx != ny) {
	  hoc_execerror("Scratch vector must be same size", 0);
	}
	for (i=0; i<nx; i++) { y[i]=x[i]; } 
        next = 1; 
	for (i=0; i < flip;) { 
	  ii = (int) ((nx+1)*drand48());
	  if (x[ii]==y[ii] && y[ii]==next) { 
	    next=x[ii]=((x[ii]==1.)?BVBASE:1.);
            i++;
	  }
	}
	return flip;
}
ENDVERBATIM
 


VERBATIM
static double vpr (void* vv) {
  int i, nx, flag, min,max;
  double* x; char c;
  FILE* f;
  nx = vector_instance_px(vv, &x);
  min=0; max=nx;
  if (ifarg(3)) {min=(int)*getarg(2); max=(int)*getarg(3)+1;} else if (ifarg(2)) {
    max=(int)*getarg(2)+1; } 
  if (min<0 || max>nx) {printf("stats vpr ERRA OOB: %d %d\n",min,max); hxe();}
  if (ifarg(1)) { 
    if (hoc_is_double_arg(1)) {
      flag=(int) *getarg(1);
      if (flag==2) { 
        for (i=min; i<max; i++) printf("%d",(x[i]>BVBASE)?1:0);
      } else if (flag==0) {
        for (i=min; i<max; i++) printf("%s%d%s",x[i]>=10?"|":"",(int)x[i],x[i]>=10?"|":"");
      } else if (flag==-1) { 
        for (i=min; i<max; i++) printf("%04.2f|",x[i]);
      } else {
        for (i=min; i<max; i++) vprpr(x[i],flag);
      }
      if (!ifarg(2)) printf("\n"); else printf(" ");
    } else {
      f = hoc_obj_file_arg(1);
      for (i=min; i<max; i++) {
        if (x[i]>BVBASE) { fprintf(f,"%d",1); 
        } else { fprintf(f,"%d",0); }
      }
      fprintf(f,"\n");
    }
  } else {
    for (i=min; i<max; i++) printf("%d",(x[i]>BVBASE)?1:0);
    printf("\n");
  }
  return 1.;
}
ENDVERBATIM



VERBATIM
static double vpr2 (void* vv) {
  int i, flag, n, nx, ny, colc, base, min,max, fl2;
  double *x, *y, cnt, ign; char c;
  nx = vector_instance_px(vv, &x);
  cnt=min=0; max=nx;
  ny = vector_arg_px(1, &y);
  if (nx!=ny) {printf("vpr2 diff sizes %d %d\n",nx,ny); hxe();}
  base=(int)*getarg(2);
  flag=(ifarg(3)?(int)*getarg(3):2);
  for (i=min,fl2=0,colc=0; i<max; i++) {
    if (flag>0 && x[i]==BVBASE && y[i]==BVBASE) {
      if (!fl2) {printf(" _ "); colc+=3;}
      fl2=1; 
    } else { 
      fl2=0; colc++;
      vprpr(x[i],base);
      if (colc>(int)newline){printf("\n    "); colc=0;}
    }
  }
  printf("\n");
  for (i=min,fl2=0,colc=0; i<max; i++) {
    if (flag>0 && x[i]==BVBASE && y[i]==BVBASE) {
      if (!fl2) {printf(" _ "); colc+=3;} 
      fl2=1; 
    } else { fl2=0;
      vprpr(y[i],base);
      colc++;
      if (colc>(int)newline){printf("\n    ",colc); colc=0;}
    }
  }
  printf("\n");
  if (flag==1) { 
    for (i=min,n=0,fl2=0,colc=0; i<max; i++) {
      if (x[i]==BVBASE && y[i]==BVBASE) { 
        fl2=1; n++;
      } else { 
        if (fl2) {printf(" %2d ",n); colc+=3;} else {printf(" "); colc++;}
        n=fl2=0; 
      }
      if (colc>(int)newline){printf("\n    ",colc); colc=0;}
    }
    printf("\n");
  }
}

static void vprpr (double x, int base) {
  int xx;
  xx=(int)x;
  if (base==0)  {    printf("%5.2f",x);
  } else if (xx>=base && base!=0) { printf("+");
  } else if (base==64) { 
    if (xx<16) {  printf("%x",xx);    
    } else if (xx<36) {  printf("%c",xx+87); 
    } else if (xx<62) {  printf("%c",xx+29); 
    } else if (xx<63) { printf("@");                
    } else if (xx<64) { printf("=");                
    } else printf("ERROR");
  } else if (base==10) {    printf("%d",xx);
  } else if (base==16) {    printf("%x",xx);
  }
}
ENDVERBATIM





VERBATIM
static double bin (void* vv) {	
  int i, j, nx, maxsz, lfl;
  double* x, *y, *ix, invl, min, max, maxf, jj;
  Object* ob;
  void* voi[2];

  min=0; max=1e9; maxf=-1e9;
  nx = vector_instance_px(vv, &x);
  ob =   *hoc_objgetarg(1);
  if (strncmp(hoc_object_name(ob),"Vector",6)==0) { lfl=0; 
    if ((maxsz=openvec(1, &y))==0) hxe();
    voi[0]=vector_arg(1);
  } else { 
    lfl=1;
    maxsz = list_vector_px3(ob, 0, &y, &voi[0]);
    if (maxsz!=(i=list_vector_px3(ob,1,&ix,&voi[1]))){printf("binERRA %d %d\n",maxsz,i); hxe();}
  }
  invl = *getarg(2);
  if (ifarg(4)) {min=*getarg(3); max=*getarg(4);
  } else if (ifarg(3)) max=*getarg(3);
  for (j=0; j<maxsz; j++) y[j]=0.;
  for (i=0; i<nx; i++) {
    if (x[i]>=max) {        y[(j=(int)(jj=(max-min)/invl))]++;
    } else if (x[i]<=min) { y[(j=0)]++; jj=0.;
    } else {
      if (x[i]>maxf) maxf=x[i]; 
      j=(int)(jj=(x[i]-min)/invl);
      if (j>maxsz-1) {printf("bin ERRB OOB: %d>%d\n",j,maxsz-1); hxe();}
      y[j]++;
    }
    if (lfl) ix[j]=jj+min;
  }
  maxsz=(max==1e9)?(int)(maxf/invl+1):(int)((max-min)/invl+1);
  vector_resize(voi[0], maxsz);
  if (lfl) vector_resize(voi[1], maxsz);
  return (double)maxsz;
}
ENDVERBATIM






VERBATIM
static double ihist (void* vv) {
  unsigned int i, j, k, n, nx, c;
  double *x, *tv, ioff, min, max, nbin, binsz; 
  ListVec* pL; Object* obl;
  nx = vector_instance_px(vv, &x); 
  i = vector_arg_px(1, &tv); 
  if (i!=nx) {printf("vecst:ihist()ERR0: diff size %d %d\n",nx,i); hxe();}
  if (!flag && !ismono1(tv,nx,1)){
    printf("vecst:ihist()ERR0A: set flag_stats for non-monotonic time vec\n",nx,i); hxe();}
  pL = AllocListVec(obl=*hoc_objgetarg(2));
  min=*getarg(3); max=*getarg(4); binsz=*getarg(5); 
  if (binsz<=0) {printf("stats:ihist()ERR0B: binsz must be >0 (%g)\n",binsz); hxe();}
  nbin=floor((max-min)/binsz);
  max=min+binsz*nbin; 
  if (verbose) printf("%g-%g in %g bins of %g\n",min,max,nbin,binsz);
  if (ifarg(6)) ioff=*getarg(6); else ioff=0.;
  if (pL->isz<2) {printf("stats:ihist()ERRA: %d\n",pL->isz); FreeListVec(&pL); hxe();}
  c=pL->isz;     
  
  for (i=0;i<c;i++) {
    pL->pv[i]=list_vector_resize(obl, i, (int)nbin);
    for (j=0;j<(int)nbin;j++) pL->pv[i][j]=0.;
  }
  i=n=0;
  if (!flag) for (;i<nx; i++) if (tv[i]>=min) break; 
  for (;i<nx; i++) { 
    if (flag) { 
      if (tv[i]>=max || tv[i]<min) continue;
    } else if (tv[i]>=max) break;
    k=(int)(x[i]-ioff); 
    if (k>=c || k<0) continue; 
    j=(int)((tv[i]-min)/binsz); 
    
    pL->pv[k][j]++;
    n++;
  }
  flag=0;
  FreeListVec(&pL);
  return (double)n;
}
ENDVERBATIM



VERBATIM
static double irate (void* vv) {
  unsigned int i, j, n, nx;
  double *prate,*phist,binsz,t1,t2;
  nx = vector_arg_px(1, &phist);
  vector_resize(vv,nx);
  vector_instance_px(vv, &prate);
  binsz = *getarg(2);
  for(i=0;i<nx;i++) {
    prate[i]=0.;
    if(phist[i]>0) {
      t1=t2=i;
      prate[i]=phist[i]*1e3/binsz;
      i++;
      break;
    }
  }
  for(;i<nx;i++) {
    prate[i]=0.;
    if(phist[i]>0) {
      t2=i;
      prate[i]=phist[i]*1e3/binsz;
      break;
    }
  }
  if(verbose>1) if(t1==t2) printf("t1==t2!\n");
  for(i=t1;i<t2;i++) prate[i] = 1e3 / (binsz*(t2-t1));
  i++;
  for(;i<nx;i++) {
    if(phist[i]>0) { 
      prate[i] = phist[i]*1e3/binsz;
      t1=t2; t2=i;      
      if(verbose>2) printf("t1 %g t2 %g\n",t1,t2);
    } else {
      if(verbose>1) if(t1==t2) printf("t1==t2!\n");
      prate[i] = 1e3 / (binsz*(t2-t1));
    }
  }
  return 1.0;
}




double dlvar (double* p, int sz) {
  int i; double s,n,d;
  s=0.0;
  for(i=0;i<sz-1;i++) {
    n = p[i]-p[i+1];
    n = n*n;
    d = p[i]+p[i+1];
    d = d*d;
    s += 3*n/d;
  }
  return s / (double) (sz-1.0);
}

static double lvar (void* vv) {
  double *x; int sz;
  if((sz=vector_instance_px(vv,&x))<2) {
    fprintf(stderr,"lvar WARNA: vector size must be >= 2, returning 0!\n");
    return 0.0;
  }
  return dlvar(x,sz);
}
ENDVERBATIM

VERBATIM





int choldc(double **a, int n, double p[])
{
  int i,j,k;
  double sum;  
  for (i=0;i<n;i++) {
    for (j=i;j<n;j++) {
      for (sum=a[i][j],k=i-1;k>=0;k--) sum -= a[i][k]*a[j][k]; 
      if (i == j) {
        if (sum <= 0.0) return 0;
        p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }
  return 1;
}
ENDVERBATIM





FUNCTION mktscor () {
VERBATIM
  double dret,*ptmp,**cF,**Y,**X,**cFT,r;
  ListVec *pTS;
  int i,j,k,nrows,ncols,nrcf;
  dret=0.0;
  pTS=0x0; ptmp=0x0; cF=cFT=0x0; X=Y=0x0;
  r = *getarg(1);
  pTS = AllocListVec(*hoc_objgetarg(2));
  nrows=pTS->plen[0]; ncols=pTS->isz;
  X = getdouble2D(nrows,ncols);
  for(i=0;i<ncols;i++) for(j=0;j<nrows;j++) X[j][i] = pTS->pv[i][j];
  nrcf=ncols;
  cF = getdouble2D(nrcf,nrcf);
  for(i=0;i<ncols;i++) for(j=0;j<=i;j++) if(j==i) cF[i][j]=1; else cF[i][j]=cF[j][i]=r;
  ptmp = (double*)calloc(ncols,sizeof(double));
  if(!choldc(cF,ncols,ptmp)) { printf("mktscor ERRA: arg must be positive definite!\n"); goto MKTSCORFREE; }
  for(i=0;i<ncols;i++) for(j=1;j<ncols;j++) if(j>i) cF[i][j]=0.0; else if(j==i) cF[i][j]=ptmp[i];
  cFT = getdouble2D(nrcf,nrcf);
  for(i=0;i<nrcf;i++) for(j=0;j<nrcf;j++) cFT[i][j] = cF[j][i];
  if(verbose){
    printf("\n\ncholsky decomp:\n");
    for(i=0;i<nrcf;i++) { for(j=0;j<nrcf;j++) printf("%g ",cFT[i][j]); printf("\n"); } printf("\n");
  }
  Y = getdouble2D(nrows,ncols);   
  for(i=0;i<nrows;i++) for(j=0;j<ncols;j++) for(k=0;k<nrcf;k++) Y[i][j] += X[i][k] * cFT[k][j];
  for(i=0;i<nrows;i++) for(j=0;j<ncols;j++) pTS->pv[j][i]=Y[i][j];
MKTSCORFREE:
  dret=1.0;
  if(pTS) FreeListVec(&pTS);
  if(ptmp) free(ptmp);
  if(cF) freedouble2D(&cF,nrcf);
  if(cFT) freedouble2D(&cFT,nrcf);
  if(Y) freedouble2D(&Y,nrows);
  if(X) freedouble2D(&X,nrows);
  return dret;
ENDVERBATIM
}


PROCEDURE install () {
  if (INSTALLED==1) {
    printf("$Id
  } else {
  INSTALLED=1
VERBATIM
  x1x=y1y=z1z=0x0;
  ihigh=1;
  install_vector_method("slope", slope);
  install_vector_method("moment", moment);
  install_vector_method("vslope", vslope);
  install_vector_method("stats", stats);
  install_vector_method("pcorrel", pcorrel);
  install_vector_method("scorrel",scorrel);
  install_vector_method("kcorrel",kcorrel);
  install_vector_method("rms",rms);
  install_vector_method("vstats", vstats);
  install_vector_method("randwd", randwd);
  install_vector_method("hamming", hamming);
  install_vector_method("flipbits", flipbits);
  install_vector_method("flipbalbits", flipbalbits);
  install_vector_method("vpr", vpr);
  install_vector_method("vpr2", vpr2);
  install_vector_method("bin", bin);
  install_vector_method("ihist", ihist);
  install_vector_method("setrnd", setrnd);
  install_vector_method("rantran", rantran);
  install_vector_method("distance", distance);
  install_vector_method("ndprd", ndprd);
  install_vector_method("hash", hash);
  install_vector_method("smash", smash);
  install_vector_method("smash1", smash1);
  install_vector_method("dpro", dpro);
  install_vector_method("unnan", unnan);
  install_vector_method("combi", combi);
  install_vector_method("shuffle", shuffle);
  install_vector_method("comb", comb);
  install_vector_method("combid", combid);
  install_vector_method("rsampsig",rsampsig);
  install_vector_method("irate",irate);
  install_vector_method("cumsum",cumsum);
  install_vector_method("lvar",lvar);
ENDVERBATIM
  }
}

PROCEDURE prhash (x) {
VERBATIM {
  long unsigned int xx;
  xx=*(long unsigned int*)(&_lx);
  printf("%16lx\n",xx);
}
ENDVERBATIM
}



FUNCTION fac (n) {
VERBATIM {
    static int ntop=4;
    static double a[101]={1.,1.,2.,6.,24.};
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
    int j,n;
    n = (int)_ln;
    if (n<0) { hoc_execerror("No negative numbers ", 0); }
    if (n>100) { 
      double x,tmp,ser;
      x = _ln;
      tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp);
      ser=1.0;
      for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
      }
      return exp(-tmp+log(2.50662827465*ser));
    } else {
      while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
      }
    return a[n];
    }
}
ENDVERBATIM
}
 


FUNCTION logfac (n) {
VERBATIM {
    static int ntop=4;
    static double a[101]={1.,1.,2.,6.,24.};
    static double cof[6]={76.18009173,-86.50532033,24.01409822,
      -1.231739516,0.120858003e-2,-0.536382e-5};
    int j,n;
    n = (int)_ln;
    if (n<0) { hoc_execerror("No negative numbers ", 0); }
    if (n>100) { 
      double x,tmp,ser;
      x = _ln;
      tmp=x+5.5;
      tmp -= (x+0.5)*log(tmp);
      ser=1.0;
      for (j=0;j<=5;j++) {
        x += 1.0;
        ser += cof[j]/x;
      }
      return (-tmp+log(2.50662827465*ser));
    } else {
      while (ntop<n) {
        j=ntop++;
        a[ntop]=a[j]*ntop;
      }
    return log(a[n]);
    }
}
ENDVERBATIM
}

VERBATIM
unsigned int hashseed2 (int na, double* x) {
  int i; union dblint xx; double y, big;
  big=pow(DBL_MAX,1./(double)(na+1)); 
  ihigh=1;
  for (i=0;i<na;i++) {
    if (x[i]==0.0) continue;
    xx.d=x[i];
    if (xx.i[0]==0) { xx.i[0]=xx.i[1]; xx.i[0]<<=4; } 
    if (xx.i[1]==0) { xx.i[1]=xx.i[0]; xx.i[1]<<=4; } 
    xx.i[0]+=(i+1); xx.i[1]+=(i+1); 
    mcell_ran4_init(&xx.i[1]);
    mcell_ran4(&xx.i[0], &y, 1, big); 
    while (y>UINT_MAX) y/=1e9; 
    ihigh*=(unsigned int)y;  
  }
  return ihigh; 
}
ENDVERBATIM

FUNCTION hashseed () {
  VERBATIM
  int i,na,nb,sf=0; double *x=0;
  if (hoc_is_double_arg(2)) {
    for (na=1;ifarg(na);na++); 
    na--;
    if (na==1) nb=2; else nb=na;
    x = (double *) malloc(sizeof(double)*nb);
    for (i=1;i<=na;i++) x[i-1]=*getarg(i);
    if (na==1) x[1]=x[0]*x[0]+13; 
    sf=1;
  } else {
    nb=(int)*getarg(1);
    x=hoc_pgetarg(2);
  } 
  hashseed2(nb,x);
  if(sf) free(x); 
  _lhashseed=(double)ihigh;
  ENDVERBATIM
}


FUNCTION vseed () {
  VERBATIM
#ifdef WIN32
  if (ifarg(1)) seed=*getarg(1); else {
    printf("TIME ACCESS NOT PRESENT IN WINDOWS\n");
    hxe();
  }
  srand48((unsigned)seed);
  set_seed(seed);
  return seed;
#else
  struct  timeval tp;
  struct  timezone tzp;
  if (ifarg(1)) seed=*getarg(1); else {
    gettimeofday(&tp,&tzp);
    seed=tp.tv_usec;
  }
  srand48((unsigned)seed);
  set_seed(seed);
  srandom(seed);
  ihigh=(unsigned int)seed;
  return seed;
#endif
  ENDVERBATIM
}



FUNCTION mc4seed () {
  VERBATIM
  double x; u_int32_t idx; unsigned int y;
  if (!ifarg(1)) {
    printf("low:%d high:%d\n",(unsigned int)ilow,(unsigned int)ihigh);
    return (double)ihigh;
  } else {
    
    x = *getarg(1);
    if (x>dmaxuint) {y=(unsigned int)x%4294967296; printf("Warning truncating %20.0f to %d (ilow)\n",x,y);x=y;}
    ilow = idx = (u_int32_t)x;
    mcell_ran4_init(idx);
  }
  if (ifarg(2)) {
    x = *getarg(2);
    if (x>dmaxuint) {y=(unsigned int)x%4294967296; printf("Warning truncating %20.0f to %d (ihigh)\n",x,y);x=y;}
    ihigh=(u_int32_t) x;
  }
  return (double)ihigh;
  ENDVERBATIM
}



FUNCTION gammln (xx) {
  VERBATIM {
    double x,tmp,ser;
    static double cof[6]={76.18009173,-86.50532033,24.01409822,-1.231739516,0.120858003e-2,-0.536382e-5};
    int j;
    x=_lxx-1.0;
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
}

FUNCTION betai(a,b,x) {
VERBATIM {
  double bt;
  double gammln(),betacf();

  if (_lx < 0.0 || _lx > 1.0) {printf("Bad x in routine BETAI\n"); hxe();}
  if (_lx == 0.0 || _lx == 1.0) bt=0.0;
  else
  bt=exp(gammln(_la+_lb)-gammln(_la)-gammln(_lb)+_la*log(_lx)+_lb*log(1.0-_lx));
  if (_lx < (_la+1.0)/(_la+_lb+2.0))
  return bt*betacf(_la,_lb,_lx)/_la;
  else
  return 1.0-bt*betacf(_lb,_la,1.0-_lx)/_lb;
 }
ENDVERBATIM
}

VERBATIM
#define ITMAX 100
#define EPS 3.0e-7
ENDVERBATIM

FUNCTION betacf(a,b,x) {
VERBATIM {
  double qap,qam,qab,em,tem,d;
  double bz,bm=1.0,bp,bpp;
  double az=1.0,am=1.0,ap,app,aold;
  int m;
  void nrerror();

  qab=_la+_lb;
  qap=_la+1.0;
  qam=_la-1.0;
  bz=1.0-qab*_lx/qap;
  for (m=1;m<=ITMAX;m++) {
    em=(double) m;
    tem=em+em;
    d=em*(_lb-em)*_lx/((qam+tem)*(_la+tem));
    ap=az+d*am;
    bp=bz+d*bm;
    d = -(_la+em)*(qab+em)*_lx/((qap+tem)*(_la+tem));
    app=ap+d*az;
    bpp=bp+d*bz;
    aold=az;
    am=ap/bpp;
    bm=bp/bpp;
    az=app/bpp;
    bz=1.0;
    if (fabs(az-aold) < (EPS*fabs(az))) return az;
  }
  printf("a or b too big, or ITMAX too small in BETACF"); return -1.;
}
ENDVERBATIM
}

FUNCTION symval() {
VERBATIM {
  Symbol *sym;
  sym = hoc_get_symbol(* hoc_pgargstr(1));
  
  return *(hoc_objectdata[sym->u.oboff]._pval); 
 }
ENDVERBATIM
}


FUNCTION tstat() {
  VERBATIM
  double r = fabs(*getarg(1));
  double N = *getarg(2);
  if(N < 2) { printf("tstat ERRA: N must be > 2!\n"); return -1; }
  return r * sqrt(N-2.)/sqrt(1.0-(r*r));
  ENDVERBATIM
}

FUNCTION tdistrib() {
  VERBATIM
  double gammln();
  double x = *getarg(1);
  double dof = *getarg(2);
  double res = (gammln( (dof+1.0) / 2.0 )  / gammln( dof / 2.0 ) );
  double pi = 3.14159265358979323846;
  res *= (1.0 / sqrt( dof * pi ) );
  res *= pow((1 + x*x/dof),-1.0*((dof+1.0)/2.0));
  return res;
  ENDVERBATIM
}









FUNCTION rcrit () {
  VERBATIM
  double rtbl[68][3]={
    1 , 0.997 , 0.999 ,    2 , 0.950 , 0.990 ,    3 , 0.878 , 0.959 ,    4 , 0.811 , 0.917 ,    5 , 0.754 , 0.874 ,
    6 , 0.707 , 0.834 ,    7 , 0.666 , 0.798 ,    8 , 0.632 , 0.765 ,    9 , 0.602 , 0.735 ,    10 , 0.576 , 0.708 ,
    11 , 0.553 , 0.684 ,    12 , 0.532 , 0.661 ,    13 , 0.514 , 0.641 ,    15 , 0.482 , 0.606 ,    16 , 0.468 , 0.590 ,
    17 , 0.456 , 0.575 ,    18 , 0.444 , 0.561 ,    19 , 0.433 , 0.549 ,    20 , 0.423 , 0.537 ,    21 , 0.413 , 0.526 ,
    22 , 0.404 , 0.515 ,    23 , 0.396 , 0.505 ,    24 , 0.388 , 0.496 ,    25 , 0.331 , 0.487 ,    26 , 0.374 , 0.478 ,
    27 , 0.367 , 0.470 ,    28 , 0.361 , 0.463 ,    29 , 0.355 , 0.456 ,    30 , 0.349 , 0.449 ,    31 , 0.344 , 0.442 ,
    32 , 0.339 , 0.436 ,    33 , 0.334 , 0.430 ,    34 , 0.329 , 0.424 ,    35 , 0.325 , 0.418 ,    36 , 0.320 , 0.413 ,
    37 , 0.316 , 0.408 ,    38 , 0.312 , 0.403 ,    39 , 0.308 , 0.398 ,    40 , 0.304 , 0.393 ,    41 , 0.301 , 0.389 ,
    42 , 0.297 , 0.384 ,    43 , 0.294 , 0.380 ,    44 , 0.291 , 0.376 ,    45 , 0.288 , 0.372 ,    46 , 0.284 , 0.368 ,
    47 , 0.281 , 0.364 ,    48 , 0.279 , 0.361 ,    58 , 0.254 , 0.330 ,    63 , 0.244 , 0.317 ,    68 , 0.235 , 0.306 ,
    73 , 0.227 , 0.296 ,    78 , 0.220 , 0.286 ,    83 , 0.213 , 0.278 ,    88 , 0.207 , 0.270 ,    93 , 0.202 , 0.263 ,
    98 , 0.195 , 0.256 ,    123 , 0.170 , 0.230 ,    148 , 0.159 , 0.210 ,    173 , 0.148 , 0.194 , 
    198 , 0.138 , 0.181 ,   298 , 0.113 , 0.148 ,    398 , 0.098 , 0.128 ,    498 , 0.088 , 0.115 ,    
    598 , 0.080 , 0.105 ,   698 , 0.074 , 0.097 ,    798 , 0.070 , 0.091 ,    898 , 0.065 , 0.086 ,  
    998 , 0.062 , 0.081   
  };
  double N;
  int get99, i , tablen, df;
  N = *getarg(1); 
  get99 = ifarg(2) ? (int)*getarg(2) : 0;
  tablen = 68;

  if(N < 3){
    printf("rcrit ERRA: N must be >= 3!\n");
    return -1.0;
  }

  if( N > 998+2 ) { printf("rcrit WARNA: Using N=1000 as estimate.\n"); N = 998+2; }

  df = (int)N - 2;

  for(i=0;i<tablen;i++) if(rtbl[i][0]==df) if(get99) return rtbl[i][2]; else return rtbl[i][1];

  for(i=1;i<tablen;i++) {
    if (rtbl[i][0] > df) {
      if(get99)
        return rtbl[i-1][2] + ((rtbl[i][2] - rtbl[i-1][2])*((df - rtbl[i-1][0] )/(rtbl[i][0] - rtbl[i-1][0])));
      else 
        return rtbl[i-1][1] + ((rtbl[i][1] - rtbl[i-1][1])*((df - rtbl[i-1][0] )/(rtbl[i][0] - rtbl[i-1][0])));      
    }
  }
  return -1.0;
  ENDVERBATIM
}