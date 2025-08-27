NEURON {
  SUFFIX sampen
  GLOBAL INSTALLED,verbose
}

PARAMETER {
  INSTALLED=0
  verbose=0
}

VERBATIM

extern int vector_instance_px();
extern void* vector_arg();
extern double *vector_newsize();



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int sampen(double *y, int mm, double r, int n, double* est,long* run,long* lastrun,double* A,double* B,double* p);
int sampen2(double *y, int mm, double r, int n, double* est, double* stdev);
void normalize(double *mdata, int n);

void getmeanstd (double* p,int n,double* mean,double* std) {
  int i;
  *mean=*std=0;
  for(i=0;i<n;i++) {
    *mean+=p[i];    *std+=p[i]*p[i];
  }
  *mean /= (double) n;
  *std = *std/(double)n - mean[0]*mean[0];
  if(*std>=0.) *std=sqrt(*std); else *std=0.0;
}


void normalize(double *mdata, int n)
{
  int i;
  double mean = 0;
  double std = 0;
  getmeanstd(mdata,n,&mean,&std);
  if(std<=0) std=1.0;
  for(i=0;i<n;i++) mdata[i] = (mdata[i]-mean) / std;
}


int sampen2(double *y, int mm, double r, int n, double* est, double* stdev)
{
    double *p = NULL;
    double *v1 = NULL, *v2 = NULL, *s1 = NULL, dv;
    int *R1 = NULL, *R2 = NULL, *F2 = NULL, *F1 = NULL, *F = NULL, FF;
    int *run = NULL, *run1 = NULL;
    double *A = NULL, *B = NULL;
    double *K = NULL, *n1 = NULL, *n2 = NULL;
    int MM;
    int m, m1, i, j, nj, jj, d, d2, i1, i2, dd;
    int nm1, nm2, nm3, nm4;
    double y1;
    int zflag=0;

    mm++;
    MM = 2 * mm;

    if ((run = (int *) calloc(n, sizeof(int))) == NULL)
        return 0;
    if ((run1 = (int *) calloc(n, sizeof(int))) == NULL)
        return 0;
    if ((R1 = (int *) calloc(n * MM, sizeof(int))) == NULL)
	return 0;
    if ((R2 = (int *) calloc(n * MM, sizeof(int))) == NULL)
	return 0;
    if ((F = (int *) calloc(n * MM, sizeof(int))) == NULL)
	return 0;
    if ((F1 = (int *) calloc(n * mm, sizeof(int))) == NULL)
	return 0;
    if ((F2 = (int *) calloc(n * mm, sizeof(int))) == NULL)
	return 0;
    if ((K = (double *) calloc((mm + 1) * mm, sizeof(double))) == NULL)
	return 0;
    if ((A = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;
    if ((B = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;
    if ((p = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;
    if ((v1 = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;
    if ((v2 = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;
    if ((s1 = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;
    if ((n1 = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;
    if ((n2 = (double *) calloc(mm, sizeof(double))) == NULL)
	return 0;

    for (i = 0; i < n - 1; i++) {
	nj = n - i - 1;
	y1 = y[i];
	for (jj = 0; jj < nj; jj++) {
	    j = jj + i + 1;
	    if (((y[j] - y1) < r) && ((y1 - y[j]) < r)) {
		run[jj] = run1[jj] + 1;
		m1 = (mm < run[jj]) ? mm : run[jj];
		for (m = 0; m < m1; m++) {
		    A[m]++;
		    if (j < n - 1)
			B[m]++;
		    F1[i + m * n]++;
		    F[i + n * m]++;
		    F[j + n * m]++;
		}
	    }
	    else
		run[jj] = 0;
	}			

	for (j = 0; j < MM; j++) {
	    run1[j] = run[j];
	    R1[i + n * j] = run[j];

	}
	if (nj > MM - 1)
	    for (j = MM; j < nj; j++)
		run1[j] = run[j];
    }				

    for (i = 1; i < MM; i++)
	for (j = 0; j < i - 1; j++)
	    R2[i + n * j] = R1[i - j - 1 + n * j];
    for (i = MM; i < n; i++)
	for (j = 0; j < MM; j++)
	    R2[i + n * j] = R1[i - j - 1 + n * j];
    for (i = 0; i < n; i++)
	for (m = 0; m < mm; m++) {
	    FF = F[i + n * m];
	    F2[i + n * m] = FF - F1[i + n * m];
	    K[(mm + 1) * m] += FF * (FF - 1);
	}

    for (m = mm - 1; m > 0; m--)
	B[m] = B[m - 1];
    B[0] = (double) n *(n - 1) / 2;
    for (m = 0; m < mm; m++) {
	p[m] = (double) A[m] / B[m];
	v2[m] = p[m] * (1 - p[m]) / B[m];
    }
    dd = 1;
    for (m = 0; m < mm; m++) {
	d2 = m + 1 < mm - 1 ? m + 1 : mm - 1;
	for (d = 0; d < d2 + 1; d++) {
	    for (i1 = d + 1; i1 < n; i1++) {
		i2 = i1 - d - 1;
		nm1 = F1[i1 + n * m];
		nm3 = F1[i2 + n * m];
		nm2 = F2[i1 + n * m];
		nm4 = F2[i2 + n * m];
		for (j = 0; j < (dd - 1); j++) {
		    if (R1[i1 + n * j] >= m + 1)
			nm1--;
		    if (R2[i1 + n * j] >= m + 1)
			nm4--;
		}
		for (j = 0; j < 2 * (d + 1); j++)
		    if (R2[i1 + n * j] >= m + 1)
			nm2--;
		for (j = 0; j < (2 * d + 1); j++)
		    if (R1[i2 + n * j] >= m + 1)
			nm3--;
		K[d + 1 + (mm + 1) * m] +=
		    (double) 2 *(nm1 + nm2) * (nm3 + nm4);
	    }
	}
    }

    n1[0] = (double) n *(n - 1) * (n - 2);
    for (m = 0; m < mm - 1; m++)
	for (j = 0; j < m + 2; j++)
	    n1[m + 1] += K[j + (mm + 1) * m];
    for (m = 0; m < mm; m++) {
	for (j = 0; j < m + 1; j++)
	    n2[m] += K[j + (mm + 1) * m];
    }

    for (m = 0; m < mm; m++) {
	v1[m] = v2[m];
	dv = (n2[m] - n1[m] * p[m] * p[m]) / (B[m] * B[m]);
	if (dv > 0)
	    v1[m] += dv;
	s1[m] = (double) sqrt((double) (v1[m]));
    }

    for (m = 0; m < mm; m++) {
	if (p[m] == 0){
            zflag=1;  
	    if(verbose>0)printf("No matches! SampEn((%d,%g,%d) = Inf"
		   " (standard deviation = Inf)!\n", m, r, n);
        } else{

            
            *est = -log(p[m]); 
            *stdev = s1[m]; 

	    if(verbose>1) printf("SampEn(%d,%g,%d) = %lf (standard deviation = %lf)\n",
		   m, r, n, *est, s1[m]);
        }
    }

    free(A);
    free(B);
    free(p);
    free(run);
    free(run1);
    free(s1);
    free(K);
    free(n1);
    free(R1);
    free(R2);
    free(v1);
    free(v2);
    free(F);
    free(F1);
    free(F2);

    if(zflag) return 2;
    return 1;
}
































int sampen (double *y, int M, double r, int n, double* est, long* run, long* lastrun, double* A, double* B, double* p)
{
    long N;
    int M1, j, nj, jj, m;
    int i;
    double y1;
    int zflag = 0;

    M++;
    memset(run,0,n*sizeof(long));
    memset(lastrun,0,n*sizeof(long));
    memset(A,0,M*sizeof(double));
    memset(B,0,M*sizeof(double));
    memset(p,0,M*sizeof(double));

    
    for (i = 0; i < n - 1; i++) {
	nj = n - i - 1;
	y1 = y[i];
	for (jj = 0; jj < nj; jj++) {
	    j = jj + i + 1;
	    if (((y[j] - y1) < r) && ((y1 - y[j]) < r)) {
		run[jj] = lastrun[jj] + 1;
		M1 = M < run[jj] ? M : run[jj];
		for (m = 0; m < M1; m++) {
		    A[m]++;
		    if (j < n - 1)
			B[m]++;
		}
	    }
	    else
		run[jj] = 0;
	}			
	for (j = 0; j < nj; j++)
	    lastrun[j] = run[j];
    }				

    N = (long) (n * (n - 1) / 2);
    p[0] = A[0] / N;

    m = M-1;
    p[m] = A[m] / B[m - 1];
    if (p[m] == 0) {
      zflag = 1;
      if(verbose>0) printf("No matches! SampEn((%d,%g,%d) = Inf!\n", m, r, n);
    } else { 
      *est = -log(p[m]);
      if(verbose>1) printf("SampEn(%d,%g,%d) = %lf\n", m, r, n, est);
    }

    if(zflag == 1) return 2; 
    return 1;
}

double* getcopy(double* in,int sz){
  double* out;
  out=(double*)malloc(sizeof(double)*sz);
  memcpy(out,in,sz*sizeof(double));
  return out;
}



static double vsampen (void* vv) {
  int n, good = 0 , getstdev = 0, sampenM;
  double* x , *outv, sampenR, sampenN, est = 0.0, stdev = 0.0, mean, std;
  long *run,*lastrun;
  double *A,*B,*p;
  if((n=vector_instance_px(vv,&x))==0){
    printf("vsampen ERRA: size 0 vector!!\n");
    return -1.0;
  }
  sampenM = ifarg(1) ? (int)*getarg(1) : 2; 
  sampenR = ifarg(2) ? *getarg(2) : 0.2; 
  sampenN = ifarg(3) ? *getarg(3) : 0.0; 
  getstdev= ifarg(4) ? *getarg(4) : 0;
  outv=ifarg(5)?vector_newsize(vector_arg(5), getstdev ? 2 : 1):0x0;  
  if(sampenN){ 
    x=getcopy(x,n); normalize(x,n);
  } else {
    mean=std=0;
    getmeanstd(x,n,&mean,&std); 
    if(std>0) sampenR = sampenR * std;
  }
  if(!getstdev) {
    run=(long*)malloc(sizeof(long)*n);   
    lastrun=(long*)malloc(sizeof(long)*n);
    A=(double*)malloc(sizeof(double)*(sampenM+1));
    B=(double*)malloc(sizeof(double)*(sampenM+1));
    p=(double*)malloc(sizeof(double)*(sampenM+1));
    good=sampen(x,sampenM,sampenR,n,&est,run,lastrun,A,B,p);
    free(run); free(lastrun); free(A); free(B); free(p);
  } else good=sampen2(x,sampenM,sampenR,n,&est,&stdev);
  if(good==2) stdev=est=-1;
  if(outv){
    outv[0]=est; if(getstdev) outv[1]=stdev;
  }
  if(sampenN) free(x);
  if(good==0){
    printf("vsampen ERRC: couldn't compute sample entropy!\n");
    return -1.0;
  }
  return est;
}





static double vsampenvst (void* vv) {
  int n, good = 0, winsz,i , j, sidx, eidx , osz , nsz, sampenM;
  double *x , *outv, sampenR, sampenN, est = 0.0, stdev = 0.0 ,*xn =0x0, mean, std;
  long *run,*lastrun;
  double *A,*B,*p;
  if((n=vector_instance_px(vv,&x))==0){
    printf("vsampenvst ERRA: size 0 vector!!\n");
    return -1.0;
  }
  sampenM = (int)*getarg(1); 
  sampenR = *getarg(2); 
  sampenN = *getarg(3); 
  if((winsz=(int)*getarg(4))<1) {  
    printf("vsampenvst ERRB: invalid window size %d!\n",winsz);
    return -1.0;
  }
  osz=ceil((double)n/winsz); if(verbose) printf("osz=%d\n",osz);
  outv=vector_newsize(vector_arg(5), osz);  

  run=(long*)malloc(sizeof(long)*winsz);   
  lastrun=(long*)malloc(sizeof(long)*winsz);
  A=(double*)malloc(sizeof(double)*(sampenM+1));
  B=(double*)malloc(sizeof(double)*(sampenM+1));
  p=(double*)malloc(sizeof(double)*(sampenM+1));

  if(sampenN) { 
    xn=getcopy(x,n); normalize(xn,n); 
    for(sidx=0,i=0;sidx<n;sidx+=winsz) {
      eidx=sidx+winsz-1; if(eidx>=n) eidx=n-1; nsz=eidx-sidx+1;
      if(verbose) if(i%20==0) printf("i:%d, sidx:%d, eidx:%d\n",i,sidx,eidx);
      good=sampen(&xn[sidx],sampenM,sampenR,nsz,&est,run,lastrun,A,B,p);
      outv[i++]=good==2?-1:est;
    }
    free(xn);
  } else {
    mean=std=0;
    getmeanstd(x,n,&mean,&std); 
    if(std>0) sampenR = sampenR * std;
    for(sidx=0,i=0;sidx<n;sidx+=winsz) {
      eidx=sidx+winsz-1; if(eidx>=n) eidx=n-1;
      if(verbose) if(i%20==0) printf("i:%d, sidx:%d, eidx:%d\n",i,sidx,eidx);
      good=sampen(&x[sidx],sampenM,sampenR,eidx-sidx+1,&est,run,lastrun,A,B,p);
      outv[i++]=good==2?-1:est;
    }
  }
  free(run); free(lastrun); free(A); free(B); free(p);
  return est;
}


ENDVERBATIM

PROCEDURE install () {
  if (INSTALLED==1){
    printf("$Id
  } else {
    INSTALLED=1
VERBATIM
  install_vector_method("vsampen", vsampen);
  install_vector_method("vsampenvst", vsampenvst);
ENDVERBATIM
  }
}