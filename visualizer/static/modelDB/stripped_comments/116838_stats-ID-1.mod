NEURON {
  SUFFIX stats
  GLOBAL  INSTALLED,seed,kmeasure,verbose,self_ok_combi,hretval
}

PARAMETER {
  
  INSTALLED=0
  kmeasure=0
  verbose=0
  self_ok_combi=0
  hretval=0
}

ASSIGNED { seed }

VERBATIM
#include <float.h>
#include <pthread.h>
#include "misc.h"
#if !defined(t)
#define _pval pval
#endif
unsigned int valseed;
static double *x1x, *y1y, *z1z;

union dblint {
  int i[2];
  double d;
};

#define VRRY 50

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


static double pcorrelsmt(double *x, double* y, int n,int nth) {
  int i,nperth,idx,rc;
  double sigxy, sigx, sigy, sigx2, sigy2, ret;
  pcorst** pp;
  pthread_t* pth;
  pthread_attr_t attr;
  ret=sigxy=sigx=sigy=sigx2,sigy2=0; 
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


static double pcorrels2(double *x, double* y, int n) {
  int i;
  double sigxy, sigx, sigy, sigx2, sigy2, dn;
  sigxy=sigx=sigy=sigx2,sigy2=0; 
  dn=(double)n;
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

static double pcorrel(void* vv) {
 int i, n;
  double *x, *y;
  n = vector_instance_px(vv, &x);
  if ((i=vector_arg_px(1, &y)) != n ) {printf("pcorrelsERRA: %d %d\n",n,i); hxe();}
  if(ifarg(2)) {
    return pcorrelsmt(x,y,n,(int)*getarg(2));
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
  rank = (double*)malloc(n*sizeof(double));
  if (!rank) return NULL;
  index = (int*)malloc(n*sizeof(int));
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
  tdata1 = (double*)malloc(n*sizeof(double));
  if(!tdata1) return 0.0; 
  tdata2 = (double*)malloc(n*sizeof(double));
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
  if ((i=vector_arg_px(1, &y)) != n ) {printf("scorrelsERRA: %d %d\n",n,i); hxe();}
  return spearman(n,x,y);
}

ENDVERBATIM
 

VERBATIM
static double unnan(void *vv) {
  int i,nx,cnt; double newnan,newinf;
  union dblint xx;
  double *x;
  newnan=newinf=0;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) newinf=newnan=*getarg(1);
  if (ifarg(2)) newinf=*getarg(2);
  for (i=0,cnt=0;i<nx;i++) { 
    xx.d=x[i];
    if (xx.i[0]==0x0 && xx.i[1]==0xfff80000) {x[i]=newnan; cnt++;}
    if (xx.i[0]==0x0 && xx.i[1]==0x7ff00000) {x[i]=newinf; cnt++;}
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
      mcell_ran4_init(xx.i[1]);
      mcell_ran4((unsigned int*)&xx.i[0], &y, 1, big); 
      prod*=y;  
    }
    if (! vfl) x[i]=prod; else return prod; 
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
    if (ifarg(i)) { y=*getarg(i++); if (y) valseed=(unsigned int)y; } 
    if (max==0) { for (i=0;i<nx;i++) x[i]=0.;
    } else mcell_ran4(&valseed, x, nx, max);
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
    return (double)valseed;
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
    if (ifarg(i)) { y=*getarg(i); if (y) valseed=(unsigned int)y; }
    if (n<=1) { for (i=0;i<nx;i++) x[i]=0.0;
    } else mcell_ran4(&valseed, x, nx, (double)n);
    if (dfl==5.5)  { for (i=0;i<nx;i++) x[i]=min+floor(x[i]);
    } else if (ex) { for (i=0;i<nx;i++) x[i]=  ex[(int)x[i]];
    } else           for (i=0;i<nx;i++) x[i]=    floor(x[i]);
    return (double)valseed;
  } else if (flag==6) { 
    min=*getarg(2); max=*getarg(3); i=4; nex=0;
    if (ifarg(i+1)) {
      nex=vector_arg_px(i, &ex); 
      valseed=(unsigned int)(*getarg(i+1));
    } else if (ifarg(i)) {
      if (hoc_is_object_arg(i)) { nex=vector_arg_px(i, &ex); 
      } else { valseed=(unsigned int)(*getarg(i)); }
    } 
    
    if (nex>1) { 
      scrset(nex);
      x1x = (double *)realloc(x1x,sizeof(double)*nx*4);
      for (i=0;i<nex;i++) scr[i]=i;
      nrn_mlh_gsort(ex, (int*)scr, nex, cmpdfn);
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
      mcell_ran4(&valseed, x1x, nx1, max-min+1-nex);
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
    mcell_ran4(&valseed, vec, 1, prn);
    if (vfl) v1[nv1]=vpr[(int)vec[0]]; else v1[nv1]=prix+floor(vec[0]);
    mcell_ran4(&valseed, vec, 1, pon);    
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
      mcell_ran4(&valseed, x1x, nx1, tot);
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


void dshuffle(double* x,int nx) {
  int n,k; double temp,y[1];
  for (n=nx;n>1;) {
    mcell_ran4(&valseed, y, 1, n);
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
    mcell_ran4(&valseed, y, 1, n);
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
    mcell_ran4(&valseed, y, 1, n);
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
  vector_resize((IvocVect*)vv,kk);
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
  g0t=g1t=0x0;
  pm=pids=0x0; 
  pcombids=0x0;
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
    printf("rsampsig ERRE: vector size (%d) < nruncombs (%d)!\n",szthis,nruncombs);
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
    
    vector_resize((IvocVect*)vhso, n0); memcpy(phso,g0t,sizeof(double)*n0); 
    hoc_call_func(pHocVecFunc,0); dm0 = hretval; 
    vector_resize((IvocVect*)vhso, n1); memcpy(phso,g1t,sizeof(double)*n1); 
    hoc_call_func(pHocVecFunc,0); dm1 = hretval; 
    hoc_pushx(dm0); hoc_pushx(dm1); hoc_call_func(pHocCompFunc,2); 
    pthis[i]=onesided?hretval:fabs(hretval); 
  }
  vector_resize((IvocVect*)vv,nruncombs);
  
  vector_resize((IvocVect*)vhso,n0); memcpy(phso,g0,sizeof(double)*n0); 
  hoc_call_func(pHocVecFunc,0); dm0 = hretval; 
  vector_resize((IvocVect*)vhso,n1); memcpy(phso,g1,sizeof(double)*n1); 
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
  for (na=1;ifarg(na);na++) {} na--; 
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
        mcell_ran4(&valseed, x1x, rvn, indx);
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
          mcell_ran4(&valseed, y, 1, ixn);  
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
  for (n=nx;n>1;) {
    mcell_ran4(&valseed, y, 1, n);
    n--;
    k=(int)y[0]; 
    temp = x[n];
    x[n] = x[k];
    x[k] = temp;
  }
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
  mcell_ran4(&valseed, y, (unsigned int)ny, (double)nx); 
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
  int i, nx, flag;
  double* x; char c;
  FILE* f;
  nx = vector_instance_px(vv, &x);
  if (ifarg(1)) { 
    if (hoc_is_double_arg(1)) {
      flag=(int) *getarg(1);
      if (flag==0) {
        for (i=0; i<nx; i++) printf("%s%d%s",x[i]>=10?"|":"",(int)x[i],x[i]>=10?"|":"");
      } else if (flag==-1) { 
        for (i=0; i<nx; i++) printf("%04.2f|",x[i]);
      } else if (flag==10) { 
        for (i=0; i<nx; i++) if (x[i]>=10) printf("+"); else printf("%d",(int)x[i]);
      } else if (flag==16) { 
        for (i=0; i<nx; i++) if (x[i]>=16) printf("+"); else printf("%x",(int)x[i]);
      } else if (flag==64) { 
        for (i=0; i<nx; i++) {
          if (x[i]>=64) {         printf("+"); 
          } else if (x[i]<16) {  printf("%x",(int)x[i]);    
          } else if (x[i]<36) {  printf("%c",(int)x[i]+87); 
          } else if (x[i]<62) {  printf("%c",(int)x[i]+29); 
          } else if (x[i]<63) { printf("@");                
          } else if (x[i]<64) { printf("=");                
          } else printf("ERROR");
        }
      }
      if (!ifarg(2)) printf("\n"); else printf(" ");
    } else {
      f = hoc_obj_file_arg(1);
      for (i=0; i<nx; i++) {
        if (x[i]>BVBASE) { fprintf(f,"%d",1); 
        } else { fprintf(f,"%d",0); }
      }
      fprintf(f,"\n");
    }
  } else {
    for (i=0; i<nx; i++) {
      if (x[i]>BVBASE) { printf("%d",1); 
      } else { printf("%d",0); }
    }
    printf("\n");
  }
  return 1.;
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
  vector_resize((IvocVect*)voi[0], maxsz);
  if (lfl) vector_resize((IvocVect*)voi[1], maxsz);
  return (double)maxsz;
}
ENDVERBATIM


PROCEDURE install () {
  if (INSTALLED==1) {
    printf("$Id
  } else {
  INSTALLED=1
VERBATIM
  x1x=y1y=z1z=0x0;
  valseed=1;
  install_vector_method("slope", slope);
  install_vector_method("moment", moment);
  install_vector_method("vslope", vslope);
  install_vector_method("stats", stats);
  install_vector_method("pcorrel", pcorrel);
  install_vector_method("scorrel",scorrel);
  install_vector_method("vstats", vstats);
  install_vector_method("randwd", randwd);
  install_vector_method("hamming", hamming);
  install_vector_method("flipbits", flipbits);
  install_vector_method("flipbalbits", flipbalbits);
  install_vector_method("vpr", vpr);
  install_vector_method("bin", bin);
  install_vector_method("setrnd", setrnd);
  install_vector_method("rantran", rantran);
  install_vector_method("distance", distance);
  install_vector_method("ndprd", ndprd);
  install_vector_method("hash", hash);
  install_vector_method("unnan", unnan);
  install_vector_method("combi", combi);
  install_vector_method("shuffle", shuffle);
  install_vector_method("comb", comb);
  install_vector_method("combid", combid);
  install_vector_method("rsampsig",rsampsig);
ENDVERBATIM
  }
}

PROCEDURE prhash (x) {
VERBATIM {
  union dblint xx;
  xx.d=_lx;
  printf("%08x%08x\n",xx.i[0],xx.i[1]);
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

FUNCTION getseed () {
  VERBATIM
  seed=(double)valseed;
  return seed;
  ENDVERBATIM
}

VERBATIM
unsigned int hashseed2 (int na, double* x) {
  int i; union dblint xx; double y, big;
  big=pow(DBL_MAX,1./(double)(na+1)); 
  valseed=1;
  for (i=0;i<na;i++) {
    if (x[i]==0.0) continue;
    xx.d=x[i];
    if (xx.i[0]==0) { xx.i[0]=xx.i[1]; xx.i[0]<<=4; } 
    if (xx.i[1]==0) { xx.i[1]=xx.i[0]; xx.i[1]<<=4; } 
    xx.i[0]+=(i+1); xx.i[1]+=(i+1); 
    mcell_ran4_init(xx.i[1]);
    mcell_ran4((unsigned int*)&xx.i[0], &y, 1, big); 
    while (y>UINT_MAX) y/=1e9; 
    valseed*=(unsigned int)y;  
  }
  return valseed; 
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
  _lhashseed=(double)valseed;
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
  valseed=(unsigned int)seed;
  return seed;
#endif
  ENDVERBATIM
}



FUNCTION mc4seed () {
  VERBATIM
  int i;
  valseed=(unsigned int)(*getarg(1));
  for (i=2;ifarg(i);i++) {
    valseed*=(unsigned int)(*getarg(i));
  }
  mcell_ran4_init(valseed); 
  return valseed;
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
  double gammln(double),betacf(double,double,double);

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
  double gammln(double);
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