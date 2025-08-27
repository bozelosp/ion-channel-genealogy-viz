NEURON {
 SUFFIX tsa
 GLOBAL INSTALLED
 GLOBAL verbose
}

PARAMETER {
  INSTALLED=0
  verbose=0
}

VERBATIM
#include "misc.h"
static double *x1x, *y1y, *z1z;

extern double dfftpow(double* x,int n,double* ppow,int powlen,int* fftlen);

#define WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))       

double *nrvector(long nl, long nh)

{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+1)*sizeof(double)));
	if (!v) { nrerror("allocation failure in vector()"); return 0x0;}
	return v-nl+1;
}

void nrfree_vector(double *v, long nl, long nh)

{
  free((char*) (v+nl-1));
}

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

static double nrsqrarg;
#define SQR(a) ((nrsqrarg=(a)) == 0.0 ? 0.0 : nrsqrarg*nrsqrarg)


void nrfour1(double mdata[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(mdata[j],mdata[i]);
			SWAP(mdata[j+1],mdata[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*mdata[j]-wi*mdata[j+1];
				tempi=wr*mdata[j+1]+wi*mdata[j];
				mdata[j]=mdata[i]-tempr;
				mdata[j+1]=mdata[i+1]-tempi;
				mdata[i] += tempr;
				mdata[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void nrspctrm(double mdata[], double p[], int m, int k)
{
  
	void nrfour1(double mdata[], unsigned long nn, int isign);
	int mm,m44,m43,m4,kk,joffn,joff,j2,j,c=0;
	double w,facp,facm,*w1,*w2,sumw=0.0,den=0.0;

	mm=m+m;
	m43=(m4=mm+mm)+3;
	m44=m43+1;
	w1=nrvector(1,m4);
	w2=nrvector(1,m);
	facm=m;
	facp=1.0/m;
	for (j=1;j<=mm;j++) sumw += SQR(WINDOW(j,facm,facp));
	for (j=1;j<=m;j++) p[j]=0.0;
	for (j=1;j<=m;j++) w2[j] = mdata[c++];
	for (kk=1;kk<=k;kk++) {
		for (joff = -1;joff<=0;joff++) {
			for (j=1;j<=m;j++) w1[joff+j+j]=w2[j];
			for (j=1;j<=m;j++) w2[j] = mdata[c++];
			joffn=joff+mm;
			for (j=1;j<=m;j++) w1[joffn+j+j]=w2[j];
		}
		for (j=1;j<=mm;j++) {
			j2=j+j;
			w=WINDOW(j,facm,facp);
			w1[j2] *= w;
			w1[j2-1] *= w;
		}
		nrfour1(w1,mm,1);
		p[1] += (SQR(w1[1])+SQR(w1[2]));
		for (j=2;j<=m;j++) {
			j2=j+j;
			p[j] += (SQR(w1[j2])+SQR(w1[j2-1])
				+SQR(w1[m44-j2])+SQR(w1[m43-j2]));
		}
		den += sumw;
	}
	den *= m4;
	for (j=1;j<=m;j++) p[j] /= den;
	nrfree_vector(w2,1,m);
	nrfree_vector(w1,1,m4);
}

double mymax(double x, double y){
  return x > y ? x : y;
}



void mysvd(int m, int n, double** u, double w[], double** v, int* ierr)

{ int i, j, k, i1, k1, l1, its;
  double c,f,h,s,x,y,z;
  int l = 0;
  double g = 0.0;
  double scale = 0.0;
  double anorm = 0.0;
  double* rv1 = (double*)malloc(n*sizeof(double));
  if (!rv1)
  { *ierr = -1;
    return;
  }
  *ierr = 0;
  
  for (i = 0; i < n; i++)
  { l = i + 1;
    rv1[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    for (k = i; k < m; k++) scale += fabs(u[k][i]);
    if (scale != 0.0)
    { for (k = i; k < m; k++)
      { u[k][i] /= scale;
        s += u[k][i]*u[k][i];
      }
      f = u[i][i];
      g = (f >= 0) ? -sqrt(s) : sqrt(s);
      h = f * g - s;
      u[i][i] = f - g;
      if (i < n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = i; k < m; k++) s += u[k][i] * u[k][j];
          f = s / h;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (k = i; k < m; k++) u[k][i] *= scale;
    }
    w[i] = scale * g;
    g = 0.0;
    s = 0.0;
    scale = 0.0;
    if (i<n-1)
    { for (k = l; k < n; k++) scale += fabs(u[i][k]);
      if (scale != 0.0)
      { for (k = l; k < n; k++)
        { u[i][k] /= scale;
          s += u[i][k] * u[i][k];
        }
        f = u[i][l];
        g = (f >= 0) ? -sqrt(s) : sqrt(s);
        h = f * g - s;
        u[i][l] = f - g;
        for (k = l; k < n; k++) rv1[k] = u[i][k] / h;
        for (j = l; j < m; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[j][k] * u[i][k];
          for (k = l; k < n; k++) u[j][k] += s * rv1[k];
        }
        for (k = l; k < n; k++)  u[i][k] *= scale;
      }
    }
    anorm = mymax(anorm,fabs(w[i])+fabs(rv1[i]));
  }
  
  for (i = n-1; i>=0; i--)
  { if (i < n-1)
    { if (g != 0.0)
      { for (j = l; j < n; j++) v[j][i] = (u[i][j] / u[i][l]) / g;
        
        for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < n; k++) s += u[i][k] * v[k][j];
          for (k = l; k < n; k++) v[k][j] += s * v[k][i];
        }
      }
    }
    for (j = l; j < n; j++)
    { v[i][j] = 0.0;
      v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  
  for (i = n-1; i >= 0; i--)
  { l = i + 1;
    g = w[i];
    if (i!=n-1)
      for (j = l; j < n; j++) u[i][j] = 0.0;
    if (g!=0.0)
    { if (i!=n-1)
      { for (j = l; j < n; j++)
        { s = 0.0;
          for (k = l; k < m; k++) s += u[k][i] * u[k][j];
          
          f = (s / u[i][i]) / g;
          for (k = i; k < m; k++) u[k][j] += f * u[k][i];
        }
      }
      for (j = i; j < m; j++) u[j][i] /= g;
    }
    else
      for (j = i; j < m; j++) u[j][i] = 0.0;
    u[i][i] += 1.0;
  }
  
  for (k = n-1; k >= 0; k--)
  { k1 = k-1;
    its = 0;
    while(1)
    
    { for (l = k; l >= 0; l--)
      { l1 = l-1;
        if (fabs(rv1[l]) + anorm == anorm) break;
        
        if (fabs(w[l1]) + anorm == anorm)
        
        { c = 0.0;
          s = 1.0;
          for (i = l; i <= k; i++)
          { f = s * rv1[i];
            rv1[i] *= c;
            if (fabs(f) + anorm == anorm) break;
            g = w[i];
            h = sqrt(f*f+g*g);
            w[i] = h;
            c = g / h;
            s = -f / h;
            for (j = 0; j < m; j++)
            { y = u[j][l1];
              z = u[j][i];
              u[j][l1] = y * c + z * s;
              u[j][i] = -y * s + z * c;
            }
          }
          break;
        }
      }
      
      z = w[k];
      if (l==k) 
      { if (z < 0.0)
        
        { w[k] = -z;
          for (j = 0; j < n; j++) v[j][k] = -v[j][k];
        }
        break;
      }
      else if (its==30)
      { *ierr = k;
        break;
      }
      else
      
      { its++;
        x = w[l];
        y = w[k1];
        g = rv1[k1];
        h = rv1[k];
        f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
        g = sqrt(f*f+1.0);
        f = ((x - z) * (x + z) + h * (y / (f + (f >= 0 ? g : -g)) - h)) / x;
        
        c = 1.0;
        s = 1.0;
        for (i1 = l; i1 <= k1; i1++)
        { i = i1 + 1;
          g = rv1[i];
          y = w[i];
          h = s * g;
          g = c * g;
          z = sqrt(f*f+h*h);
          rv1[i1] = z;
          c = f / z;
          s = h / z;
          f = x * c + g * s;
          g = -x * s + g * c;
          h = y * s;
          y = y * c;
          for (j = 0; j < n; j++)
          { x = v[j][i1];
            z = v[j][i];
            v[j][i1] = x * c + z * s;
            v[j][i] = -x * s + z * c;
          }
          z = sqrt(f*f+h*h);
          w[i1] = z;
          
          if (z!=0.0)
          { c = f / z;
            s = h / z;
          }
          f = c * g + s * y;
          x = -s * g + c * y;
          for (j = 0; j < m; j++)
          { y = u[j][i1];
            z = u[j][i];
            u[j][i1] = y * c + z * s;
            u[j][i] = -y * s + z * c;
          }
        }
        rv1[l] = 0.0;
        rv1[k] = f;
        w[k] = x;
      }
    }
  }
  free(rv1);
  return;
}

int mycompare( const void* a, const void* b ) {
  double* arg1 = (double*) a;
  double* arg2 = (double*) b;
  if( *arg1 < *arg2 ) return -1;
  else if( *arg1 == *arg2 ) return 0;
  else return 1;
}  

double* myspectrum(double* v1,int dc,int* anslen){
  
  
  

  
  

  
  int mr;
  
  mr = dc / 8;

  
  int m = 1;
  while(m<mr) m*=2;
  *anslen = m;

  int k = (int) (ceil(((double)(dc)/m-1.)/2.));
  int n = (2*k+1)*m;

  int i;

  double *mdata = (double *)calloc(n,(unsigned)sizeof(double));
  for (i=0;i<dc;++i) mdata[i] = v1[i]; 

  
  double* ans = (double *)calloc(m,(unsigned)sizeof(double));
  
  nrspctrm(&mdata[0], &ans[0], m, k);

  free((char *)mdata);

  return ans;
}

double myspct(void* v){
  double* x;
  int n = vector_instance_px(v,&x) , anslen = 0 , i = 0;
  double* ans = myspectrum(x,n,&anslen);
  for(i=0;i<anslen;i++) printf("%g ",ans[i]);
  printf("\n");
  double* p;
  int outlen = vector_arg_px(1,&p);
  for(i=0;i<outlen && i<anslen;i++) p[i]=ans[i];
  free(ans);
  return 1.0;
}

int iinrange(double* x,int sz,double dmin,double dmax){
  int i = 0, cnt = 0;
  for(i=0;i<sz;i++) if(x[i]>=dmin && x[i]<=dmax) cnt++;
  return cnt;
}

double inrange(void* v){
  double* x , dmin, dmax;
  int sz = vector_instance_px(v,&x) , i , cnt = 0;
  return (double) iinrange(x,sz,*getarg(1),*getarg(2));
}

extern double dfftpow(double* x,int n,double* ppow,int powlen,int* fftlen);

ENDVERBATIM








FUNCTION tscoravgband () {
  VERBATIM

  int i;

  double dRet = 0.0;
  ListVec* pTS = 0;

  double* pcvin = 0;
  double** pfilter = 0x0;

  double** pTSN = 0; 

  double* vTCor = 0; 
  int corsz = vector_arg_px(1,&vTCor);

  double* pTMP = 0x0;

  ListVec* pTSOut = 0;

  pTS = AllocListVec(*hoc_objgetarg(2)); 
  
  if(!pTS){
    printf("tscoravg ERRA: problem initializing 1st 2 args!\n");
    goto CLEANUP;
  }
  int iChans = pTS->isz; 
  if(iChans < 2) {
    printf("tscoravg ERRB: must have at least 2 EEG channels!\n");
    goto CLEANUP;
  }
  int iWinSz = (int) *getarg(3); 
  int iSamples = pTS->plen[0]; 
  int iINC = ifarg(4)?(int)*getarg(4):1; 
  double *phz = 0, *psat = 0, *pfft = 0; int sz = 0 , iNoiseTest = 0;
  if(ifarg(5) && ifarg(6)){ 
    if((sz=vector_arg_px(5,&phz))!=corsz){
      printf("tscoravg ERRG: invalid size for phz should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else if((sz=vector_arg_px(6,&psat))!=corsz){
      printf("tscoravg ERRH: invalid size for psat should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else {	      
      iNoiseTest = 1;
      pfft = (double*)calloc(iWinSz+1,sizeof(double));
    }
  }
  int iOutSz = (iSamples - iWinSz + 1) / iINC;

  double *pdelta=0,*ptheta=0,*pbeta=0,*pgamma=0,*pripple=0;
  if(ifarg(7) && vector_arg_px(7,&pdelta)<iOutSz){ printf("bad pdelta size\n"); goto CLEANUP; } 
  if(ifarg(8) && vector_arg_px(8,&ptheta)<iOutSz){ printf("bad ptheta size\n"); goto CLEANUP; } 
  if(ifarg(9)&& vector_arg_px(9,&pbeta)<iOutSz){ printf("bad pbeta size\n"); goto CLEANUP; } 
  if(ifarg(10)&& vector_arg_px(10,&pgamma)<iOutSz){ printf("bad pgamma size\n"); goto CLEANUP; } 
  if(ifarg(11)&& vector_arg_px(11,&pripple)<iOutSz){ printf("bad pripple size\n"); goto CLEANUP; } 

  double sampr = *getarg(12);

  double *plohz=0x0,*phihz=0x0;
  int iFilters = 1;
  if((iFilters=vector_arg_px(13,&plohz))!=(i=vector_arg_px(14,&phihz))) {
    printf("invalid arg 13,14 lohz/hihz filter sizes: %d %d\n",iFilters,i); goto CLEANUP;
  } else if(iFilters<1) { printf("no filters specified!\n"); goto CLEANUP; }
  int N = 1, M = 1025;
  while(N < iWinSz) N*=2;
  if(N>M) M=N+1;
  pfilter = getdouble2D(iFilters,N);
  
  for(i=0;i<iFilters;i++){ plohz[i]/=sampr; phihz[i]/=sampr; }
  for(i=0;i<iFilters;i++){ wsfirBP(pfilter[i], M, 1, plohz[i], phihz[i], N); wrap(pfilter[i],N,M); }
  pcvin = (double*) calloc(N,sizeof(double));

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscoravg ERRC: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto CLEANUP;
    }
  }

  if(corsz < iOutSz){
    printf("tscoravg ERRD: need output vector of at least %d as arg 1, have only %d!\n",iOutSz,corsz);
    goto CLEANUP;
  }
  if(verbose) printf("iOutSz=%d\n",iOutSz);

  pTSN = getdouble2D(iChans,2*N); if(verbose) printf("got pTSN\n");
  if(!pTSN){
    printf("tscor ERRD: out of memory!\n");
    goto CLEANUP;
  }

  if(ifarg(15) && (pTSOut = AllocListVec(*hoc_objgetarg(15)))) { 
    if(pTSOut->isz!=iChans) { printf("incorrect output vec-list size: %d %d\n",pTSOut->isz,iChans); goto CLEANUP; }
    for(i=0;i<iChans;i++) if(pTSOut->plen[i] < iSamples) { printf("incorrect output vec-size %d %d\n",pTSOut->plen[i],iSamples);
      goto CLEANUP; }}

  pTMP = (double*) calloc(2*N,sizeof(double));

  double* psums = (double*)calloc(iChans,sizeof(double)); 
  double* psums2 = (double*)calloc(iChans,sizeof(double));

  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iOutIDX = 0 , iFilt = 0;
  for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) vTCor[iOutIDX] = 0.;
  if( phz ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) phz[iOutIDX]=0.;
  if( psat ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) psat[iOutIDX]=0.;
  iOutIDX=0;
  double CP = iChans*(iChans-1.0)/2.0;
  for(iStartIDX=0;iOutIDX<iOutSz && iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iOutIDX++){
    if(iOutIDX%1000==0) printf("%d\n",iOutIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    if(phz) phz[iOutIDX]=0.; if(psat)psat[iOutIDX]=0.;
    if(pdelta) pdelta[iOutIDX]=0.; if(ptheta) ptheta[iOutIDX]=0.;
    if(pbeta) pbeta[iOutIDX]=0.; if(pgamma) pgamma[iOutIDX]=0.; if(pripple) pripple[iOutIDX]=0.;

    double dsatavg = 0., dhzavg = 0.;
    for(iChan=0;iChan<iChans;iChan++){       
      i=0;      
      int iSat = 0; 
      double* pc = &pTS->pv[iChan][iStartIDX]; 
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,pc++)if((*pc>=900. && *pc<=1010.) || (*pc>=-1010. && *pc<=-900.)) iSat++;
      dsatavg += (double)iSat/ (double)iWinSz;
      int fftlen = 0; memset(pfft,0,sizeof(double)*(iWinSz+1));
      if(dfftpow(&pTS->pv[iChan][iStartIDX],iWinSz,pfft,iWinSz+1,&fftlen)){       
	int f = 0; pfft[0]=0.0; double fsum = 0.0; double* ppf = pfft;
	for(f=0;f<fftlen && f<=sampr/2;f++) fsum += *ppf++; 
	for(f=55;f<=65 && f<fftlen;f++) dhzavg += pfft[f] / fsum; 
	for(f=115;f<=125 && f<fftlen;f++) dhzavg += pfft[f] / fsum; 
	if(pdelta) {
	  fsum=0.0;
	  for(f=0;f<fftlen;f++) fsum+=pfft[f];
	  for(f=0;f<=3;f++) pdelta[iOutIDX] += pfft[f]/fsum;
	  if(ptheta) for(f=4;f<=12;f++) ptheta[iOutIDX] += pfft[f]/fsum;
	  if(pbeta) for(f=13;f<=29;f++) pbeta[iOutIDX] += pfft[f]/fsum;
	  if(pgamma) for(f=30;f<=100;f++) pgamma[iOutIDX] += pfft[f]/fsum;
	  if(pripple) for(f=101;f<=300 && f<fftlen;f++) pripple[iOutIDX] += pfft[f]/fsum;
	}
      }	
    }
    phz[iOutIDX] = dhzavg / (double)iChans;
    psat[iOutIDX] = dsatavg / (double)iChans;
    if(pdelta) pdelta[iOutIDX] /= (double)iChans;
    if(ptheta) ptheta[iOutIDX] /= (double)iChans;
    if(pbeta) pbeta[iOutIDX] /= (double)iChans;
    if(pgamma) pgamma[iOutIDX] /= (double)iChans;
    if(pripple) pripple[iOutIDX] /= (double)iChans;

    
    for(iChan=0;iChan<iChans;iChan++){
      memset(pTMP,0,sizeof(double)*2*N);
      double* pC = &pTS->pv[iChan][iStartIDX], *pIN = pcvin;
      if(0) for(IDX=iStartIDX;IDX<iEndIDX;IDX++) *pIN++ = *pC++; 
      for(iFilt=0;iFilt<iFilters;iFilt++) {        
        if(1) {
          memset(pTSN[iChan],0,sizeof(double)*N);
          pIN = pcvin; pC = &pTS->pv[iChan][iStartIDX];
          for(IDX=iStartIDX;IDX<iEndIDX;IDX++) *pIN++ = *pC++; 
        }
        convlv(pcvin-1,N,pfilter[iFilt]-1,M,1,pTSN[iChan]-1);
        for(IDX=0;IDX<N;IDX++) pTMP[IDX] += pTSN[iChan][IDX]; 
      }
      for(IDX=0;IDX<N;IDX++) pTSN[iChan][IDX] = pTMP[IDX]; 
      if(pTSOut)for(IDX=0;IDX<N && IDX+iStartIDX<pTSOut->plen[iChan];IDX++) pTSOut->pv[iChan][iStartIDX+IDX]=pTSN[iChan][IDX];

      psums[iChan]=psums2[iChan]=0.;
      pC = pTSN[iChan]; 
      for(i=0;i<N;i++,pC++){ 
	psums[iChan] += *pC;
	psums2[iChan] += (*pC * *pC); 
      }
      double davg = psums[iChan] / N;
      double dstdev = psums2[iChan]/N - davg*davg; 
      if(dstdev > 0. ) dstdev = sqrt(dstdev); else dstdev=1.; 
      if(dstdev <= 0.) dstdev = 1.0;
      pC = pTSN[iChan]; 
      for(i=0;i<N;i++,pC++) *pC = (*pC - davg)/dstdev; 
    }
    
    int iC1,iC2;
    double dpsum = 0.0;
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=0;iC2<iC1;iC2++){
        double r = 0.0, *p1 = pTSN[iC1], *p2 = pTSN[iC2];
	for(i=0;i<N;i++) r += (*p1++ * *p2++);
        r /= (double) N;
        dpsum += r;
      }
    }
    dpsum /= CP;
    if(dpsum > 1.0 || dpsum < -1.0) printf("out of bounds = %g!\n",dpsum);
    vTCor[iOutIDX] = dpsum;
  }
  dRet = 1.0;
  printf("\n");

CLEANUP:
  if(pTS) FreeListVec(&pTS);                       if(pTSN) freedouble2D(&pTSN,iChans);
  free(psums);                                     free(psums2); 
  if(pfft) free(pfft);
  if(pfilter) freedouble2D(&pfilter,iFilters);     if(pcvin) free(pcvin);
  if(pTMP) free(pTMP);                             if(pTSOut) FreeListVec(&pTSOut);
  return dRet;
  ENDVERBATIM
}



FUNCTION tsnorm () {
  VERBATIM
  double dRet = 0.0;
  ListVec* pEIG = AllocListVec(*hoc_objgetarg(1));
  if(!pEIG || pEIG->isz < 1 || pEIG->plen[0]<1){
    printf("tsnorm ERRA: problem initializing 1st arg!\n");
    goto TSNCLEANUP;
  }
  double* dmean = 0, *dstdev = 0;
  int cols = vector_arg_px(2,&dmean);
  if(cols!=vector_arg_px(3,&dstdev) || cols!=pEIG->plen[0]){
    printf("tsnorm ERRB: problem initializing 2nd,3rd arg!\n");
    goto TSNCLEANUP;
  }
  int i,j;
  double val;
  for(i=0;i<pEIG->isz;i++)for(j=0;j<cols;j++)pEIG->pv[i][j]=(pEIG->pv[i][j]-dmean[j])/dstdev[j];
  dRet = 1.0;
TSNCLEANUP:
  if(pEIG) FreeListVec(&pEIG);
  return dRet;
  ENDVERBATIM
}

VERBATIM
int ddiff(double* pin,double* pout,int sz){
  int i;
  for(i=0;i<sz-1;i++) pout[i]=pin[i+1]-pin[i];
  return 1;
}
ENDVERBATIM










FUNCTION tsgetmeans () {
VERBATIM
  int iC1,iC2,iChans,idx,jdx;
  double dRet = 0.0;
  ListVec *plvcor = 0, *plvhz = 0, *plvsat = 0, *plvderiv = 0;
  double* pavgA=0,*pavgB=0,*pavgC=0,*pprctbadch=0,*pfid=0,*pbadc=0;
  int vsz = 0, fid=-1; 
  
  
  
  
  if(!(plvcor = AllocListVec(*hoc_objgetarg(1)))){
    printf("tsgetmeans ERRA: couldn't get cor vec list!\n");
    goto MCLEANUP;
  }
  if(!(plvhz = AllocListVec(*hoc_objgetarg(2)))){
    printf("tsgetmeans ERRB: couldn't get phz vec list!\n");
    goto MCLEANUP;
  }
  if(plvhz->isz < plvcor->isz){
    printf("tsgetmeans ERRC: phz vec list size < cor vec list size : %d %d\n",plvhz->isz,plvcor->isz);
    goto MCLEANUP;
  }
  if(!(plvsat = AllocListVec(*hoc_objgetarg(3)))){
    printf("tsgetmeans ERRD: couldn't get psat vec list!\n");
    goto MCLEANUP;
  }
  if(plvsat->isz < plvcor->isz){
    printf("tsgetmeans ERRE: psat vec list size < cor vec list size : %d %d\n",plvsat->isz,plvcor->isz);
    goto MCLEANUP;
  }  
  if(!(plvderiv = AllocListVec(*hoc_objgetarg(4)))){
    printf("tsgetmeans ERRF: couldn't get pderiv vec list!\n");
    goto MCLEANUP;
  }
  if(plvderiv->isz < plvcor->isz){
    printf("tsgetmeans ERRG: pderiv vec list size < cor vec list size : %d %d\n",plvderiv->isz,plvcor->isz);
    goto MCLEANUP;
  }  
  if(vector_arg_px(5,&pavgA) < plvcor->isz || 
     vector_arg_px(6,&pavgB) < plvcor->isz || 
     vector_arg_px(7,&pavgC) < plvcor->isz ||
     vector_arg_px(13,&pfid) < plvcor->isz ||
     vector_arg_px(15,&pprctbadch)< plvcor->isz){
    printf("tsgetmeans ERRH: output vecs must have size >= %d!\n",plvcor->isz);
    goto MCLEANUP;
  }
  double phzt = *getarg(8), psatt = *getarg(9) , pderivt = *getarg(10);
  iChans = (int)*getarg(11);
  double pfctr = *getarg(12);
  fid = (int)*getarg(14);
  pbadc = (double*) malloc(iChans*sizeof(double));
  
  
  

  
  
  double dSumA,dSumB,dSumC;
  int cntA,cntB,cntC,cntBadCh;
  for(idx=0;idx<plvcor->isz;idx++){
    if(fid>=0 && pfid[idx]!=fid) continue;
    cntA=cntB=cntC=jdx=0;
    dSumA=dSumB=dSumC=0.;
    for(iC1=0;iC1<iChans;iC1++){ 
      if(plvsat->pv[idx][iC1] >= psatt ||
         plvderiv->pv[idx][iC1] >= pderivt)
        pbadc[iC1]=1;
      else
        pbadc[iC1]=0;
    }
    for(iC1=0;iC1<iChans;iC1++){
      if(pbadc[iC1])continue;
      for(iC2=0;iC2<iC1;iC2++,jdx++){
        dSumA += plvcor->pv[idx][jdx];
        cntA++;
        if(pbadc[iC2]) continue;
        dSumB += plvcor->pv[idx][jdx];
        cntB++;
      }
    }
    pavgA[idx] = dSumA / (double) cntA;
    if(cntB) pavgB[idx] = dSumB / (double) cntB; else pavgB[idx] = -666.;
    jdx=0;
    for(iC1=0;iC1<iChans;iC1++) if(plvhz->pv[idx][iC1] >= phzt) pbadc[iC1]=1; 
    for(iC1=0;iC1<iChans;iC1++){
      if(pbadc[iC1]) continue;
      for(iC2=0;iC2<iC1;iC2++,jdx++){
        if(pbadc[iC2]) continue;
        dSumC += plvcor->pv[idx][jdx];
        cntC++;
      }
    }
    if(cntC) pavgC[idx] = dSumC / (double) cntC; else pavgC[idx] = -666.;
    cntBadCh=0;
    for(iC1=0;iC1<iChans;iC1++) if(pbadc[iC1]) cntBadCh++;
    pprctbadch[idx]= (double)cntBadCh/(double)iChans;
  }

  dRet = 1.0;
MCLEANUP:
 if(pbadc) free(pbadc);
 return dRet;
ENDVERBATIM
}



FUNCTION tscorfull () {
  VERBATIM

  int i;

  double dRet = 0.0;
  ListVec* pTS = 0, *plvcor = 0, *plvhz = 0, *plvsat = 0, *plvderiv = 0;

  double** pTSN = 0; 
  double* pdiff = 0; 
  double* psums = 0, *psums2 = 0, *pfft = 0;

  int corsz=0; 
  plvcor = AllocListVec(*hoc_objgetarg(1)); 
  if(!plvcor || (corsz=plvcor->isz)<1){
    printf("tscorfull ERRA: problem initializing arg1 plvcor!\n");
    goto FCLEANUP;
  }

  pTS = AllocListVec(*hoc_objgetarg(2)); 
  if(!pTS){
    printf("tscorfull ERRA: problem initializing arg2 pTS!\n");
    goto FCLEANUP;
  }

  int iChans = pTS->isz; 
  if(iChans < 2) {
    printf("tscorfull ERRB: must have at least 2 timeseries!\n");
    goto FCLEANUP;
  }

  int iMatSz = iChans*(iChans-1)/2;

  for(i=0;i<corsz;i++){   
    if(plvcor->plen[i]<iMatSz){
      printf("tscorfull ERRC: each cor vec must have sz >= %d, [%d] has only %d\n",iMatSz,i,plvcor->plen[i]);
      goto FCLEANUP;
    }
  }

  int iWinSz = (int) *getarg(3); 
  int iINC = ifarg(4)?(int)*getarg(4):1;  

  plvhz = AllocListVec(*hoc_objgetarg(5));
  plvsat = AllocListVec(*hoc_objgetarg(6));
  plvderiv = AllocListVec(*hoc_objgetarg(7));

  pfft = (double*)calloc(iWinSz+1,sizeof(double));
  pdiff = (double*)calloc(iWinSz,sizeof(double));

  int iSamples = pTS->plen[0]; 
  int iOutSz = (iSamples - iWinSz + 1) / iINC;

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscorfull ERRD: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto FCLEANUP;
    }
  }

  if(corsz < iOutSz){
    printf("tscorfull ERRE: need output list vector of at least sz %d as arg 1, have only %d!\n",iOutSz,corsz);
    goto FCLEANUP;
  }   if(verbose) printf("iOutSz=%d\n",iOutSz);
  if(plvhz->isz<iOutSz || plvsat->isz<iOutSz || plvderiv->isz<iOutSz){
    printf("tscorfull ERRF: need output plvhz,plvsat,plvderiv size of at least %d!\n",iOutSz);
    goto FCLEANUP;
  }
  for(i=0;i<iOutSz;i++){
    if(plvhz->plen[i]<iChans || plvsat->plen[i]<iChans || plvderiv->plen[i]<iChans){
      printf("tscorfull ERRG: each vec in plvhz,plvsat,plvderiv must be sz >= %d!\n",iChans);
      goto FCLEANUP;
    }
  }

  pTSN = getdouble2D(iChans,iWinSz); if(verbose) printf("got pTSN\n"); 
  if(!pTSN){
    printf("tscor ERRD: out of memory!\n");
    goto FCLEANUP;
  }

  psums = (double*)calloc(iChans,sizeof(double));   
  psums2 = (double*)calloc(iChans,sizeof(double));  
  double* pc = 0;
  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iOutIDX = 0;
  FillListVec(plvcor,0.0); FillListVec(plvhz,0.); FillListVec(plvsat,0.); FillListVec(plvderiv,0.);
  for(iStartIDX=0;iOutIDX<iOutSz && iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iOutIDX++){
    if(iOutIDX%1000==0) printf("%d\n",iOutIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    for(iChan=0;iChan<iChans;iChan++){ 
      psums[iChan]=psums2[iChan]=0.;
      pc = &pTS->pv[iChan][iStartIDX];
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,pc++){
        psums[iChan] += *pc; 
        psums2[iChan] += (*pc * *pc);
      }
    }
    for(iChan=0;iChan<iChans;iChan++){       
      i=0;
      double davg = psums[iChan] / iCurSz;
      double dstdev = psums2[iChan]/iCurSz - davg*davg; 
      if(dstdev > 0. ) dstdev = sqrt(dstdev); else dstdev=1.; 
      if(dstdev <= 0.) dstdev = 1.0;
      
      int iSat = 0;
      pc = &pTS->pv[iChan][iStartIDX];
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++,pc++){	
        if((*pc>=900. && *pc<=1010.) || (*pc>=-1010. && *pc<=-900.)) iSat++;
        pTSN[iChan][i] = (*pc - davg) / dstdev;	
      }
      plvsat->pv[iOutIDX][iChan] = (double)iSat/ (double)iWinSz;
      
      int fftlen = 0; memset(pfft,0,sizeof(double)*(iWinSz+1));
      if(dfftpow(&pTS->pv[iChan][iStartIDX],iWinSz,pfft,iWinSz+1,&fftlen)){
        int f = 0; pfft[0]=0.0; double fsum = 0.0;
        for(f=0;f<fftlen && f<=1000;f++) fsum += pfft[f]; 
        pc = &plvhz->pv[iOutIDX][iChan];
        for(f=55;f<=65 && f<fftlen;f++)   *pc += pfft[f]; 
        for(f=115;f<=125 && f<fftlen;f++) *pc += pfft[f]; 
        *pc /= fsum;
      }
      
      ddiff(&pTS->pv[iChan][iStartIDX],pdiff,iWinSz);
      plvderiv->pv[iOutIDX][iChan] = (double)iinrange(pdiff, iWinSz-1, -1.0, 1.0) / (double)(iWinSz-1.0);
    }

    
    int iC1,iC2; pc = &plvcor->pv[iOutIDX][0];
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=iC1+1;iC2<iChans;iC2++){
	for(i=0;i<iWinSz;i++) 
          *pc += pTSN[iC1][i]*pTSN[iC2][i];
        *pc /= (double) iWinSz;
        pc++;
      }
    }
  }
  dRet = 1.0;
  printf("\n");

FCLEANUP:
  if(pTS) FreeListVec(&pTS);
  if(pTSN) freedouble2D(&pTSN,iChans);
  if(psums) free(psums);
  if(psums2) free(psums2);
  if(pfft) free(pfft);
  if(pdiff) free(pdiff);
  return dRet;
  ENDVERBATIM
}


FUNCTION tscoravg () {
  VERBATIM

  int i;

  double dRet = 0.0;
  ListVec* pTS = 0;

  double** pTSN = 0; 
  double** pCorrel = 0; 

  double* vTCor = 0; 
  int corsz = vector_arg_px(1,&vTCor);

  pTS = AllocListVec(*hoc_objgetarg(2)); 
  
  if(!pTS){
    printf("tscoravg ERRA: problem initializing 1st 2 args!\n");
    goto CLEANUP;
  }

  int iChans = pTS->isz; 
  if(iChans < 2) {
    printf("tscoravg ERRB: must have at least 2 EEG channels!\n");
    goto CLEANUP;
  }
  int iWinSz = (int) *getarg(3); 


  int iSamples = pTS->plen[0]; 

  int iINC = ifarg(4)?(int)*getarg(4):1; 

  double *phz = 0, *psat = 0, *pfft = 0; int sz = 0 , iNoiseTest = 0;
  if(ifarg(5) && ifarg(6)){ 
    if((sz=vector_arg_px(5,&phz))!=corsz){
      printf("tscoravg ERRG: invalid size for phz should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else if((sz=vector_arg_px(6,&psat))!=corsz){
      printf("tscoravg ERRH: invalid size for psat should be %d, is %d!\n",corsz,sz);
      goto CLEANUP;
    } else {	      
      iNoiseTest = 1;
      pfft = (double*)calloc(iWinSz+1,sizeof(double));
    }
  }

  int idoabs = ifarg(7)?(int)*getarg(7):0; 

  int iOutSz = (iSamples - iWinSz + 1) / iINC;

  double *pdelta=0,*ptheta=0,*pbeta=0,*pgamma=0,*pripple=0;
  if(ifarg(8) && vector_arg_px(8,&pdelta)<iOutSz){ printf("bad pdelta size\n"); goto CLEANUP; } 
  if(ifarg(9) && vector_arg_px(9,&ptheta)<iOutSz){ printf("bad ptheta size\n"); goto CLEANUP; } 
  if(ifarg(10)&& vector_arg_px(10,&pbeta)<iOutSz){ printf("bad pbeta size\n"); goto CLEANUP; } 
  if(ifarg(11)&& vector_arg_px(11,&pgamma)<iOutSz){ printf("bad pgamma size\n"); goto CLEANUP; } 
  if(ifarg(12)&& vector_arg_px(12,&pripple)<iOutSz){ printf("bad pripple size\n"); goto CLEANUP; } 

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscoravg ERRC: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto CLEANUP;
    }
  }

  if(corsz < iOutSz){
    printf("tscoravg ERRD: need output vector of at least %d as arg 1, have only %d!\n",iOutSz,corsz);
    goto CLEANUP;
  }
  if(verbose) printf("iOutSz=%d\n",iOutSz);

  pCorrel = getdouble2D(iChans,iChans); if(verbose) printf("got pCorrel\n");
  pTSN = getdouble2D(iChans,iWinSz); if(verbose) printf("got pTSN\n");
  if(!pCorrel || !pTSN){
    printf("tscor ERRD: out of memory!\n");
    goto CLEANUP;
  }

  double* psums = (double*)calloc(iChans,sizeof(double));
  double* psums2 = (double*)calloc(iChans,sizeof(double));

  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iOutIDX = 0;
  for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) vTCor[iOutIDX] = 0.;
  if( phz ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) phz[iOutIDX]=0.;
  if( psat ) for(iOutIDX=0;iOutIDX<corsz;iOutIDX++) psat[iOutIDX]=0.;
  iOutIDX=0;
  for(iStartIDX=0;iOutIDX<iOutSz && iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iOutIDX++){
    if(iOutIDX%1000==0) printf("%d\n",iOutIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    if(phz) phz[iOutIDX]=0.; if(psat)psat[iOutIDX]=0.;
    if(pdelta) pdelta[iOutIDX]=0.; if(ptheta) ptheta[iOutIDX]=0.;
    if(pbeta) pbeta[iOutIDX]=0.; if(pgamma) pgamma[iOutIDX]=0.; if(pripple) pripple[iOutIDX]=0.;
    if(iStartIDX==0 || iINC > 1){
      for(iChan=0;iChan<iChans;iChan++){ 
	psums[iChan]=psums2[iChan]=0.;
	for(IDX=iStartIDX;IDX<iEndIDX;IDX++){
	  dVal = pTS->pv[iChan][IDX];
	  psums[iChan] += dVal;
	  psums2[iChan] += dVal*dVal;
	}
      }
    } else {
      for(iChan=0;iChan<iChans;iChan++){
	dVal = pTS->pv[iChan][iStartIDX-1];
	psums[iChan] -= dVal;
	psums2[iChan] -= dVal*dVal;
	dVal = pTS->pv[iChan][iEndIDX-1];
	psums[iChan] += dVal;
	psums2[iChan] += dVal*dVal;
      }
    }
    double dsatavg = 0., dhzavg = 0.;
    for(iChan=0;iChan<iChans;iChan++){       
      i=0;
      double davg = psums[iChan] / iCurSz;
      double dstdev = psums2[iChan]/iCurSz - davg*davg; 
      if(dstdev > 0. ) dstdev = sqrt(dstdev); else dstdev=1.; 
      if(dstdev <= 0.) dstdev = 1.0;
      if(iNoiseTest){	
	
	int iSat = 0;
	for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++){	
	  dVal = pTS->pv[iChan][IDX];
	  if((dVal>=900. && dVal<=1010.) || (dVal>=-1010. && dVal<=-900.)) iSat++;
	  pTSN[iChan][i] = (dVal - davg) / dstdev;	
	}
	dsatavg += (double)iSat/ (double)iWinSz;
	
	int fftlen = 0; memset(pfft,0,sizeof(double)*(iWinSz+1));
	if(dfftpow(&pTS->pv[iChan][iStartIDX],iWinSz,pfft,iWinSz+1,&fftlen)){
	  int f = 0; pfft[0]=0.0; double fsum = 0.0;
	  for(f=0;f<fftlen && f<=1000;f++) fsum += pfft[f]; 
	  for(f=55;f<=65 && f<fftlen;f++) dhzavg += pfft[f] / fsum; 
	  for(f=115;f<=125 && f<fftlen;f++) dhzavg += pfft[f] / fsum; 
	  if(pdelta) {
	    for(f=0;f<=3;f++) pdelta[iOutIDX] += pfft[f]/fsum;
	    if(ptheta) for(f=4;f<=12;f++) ptheta[iOutIDX] += pfft[f]/fsum;
	    if(pbeta) for(f=13;f<=29;f++) pbeta[iOutIDX] += pfft[f]/fsum;
	    if(pgamma) for(f=30;f<=100;f++) pgamma[iOutIDX] += pfft[f]/fsum;
	    if(pripple) for(f=101;f<=300 && f<fftlen;f++) pripple[iOutIDX] += pfft[f]/fsum;
	  }
	}	
      } else {
	for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++){	
	  dVal = pTS->pv[iChan][IDX];
	  pTSN[iChan][i] = (dVal - davg) / dstdev;	
	}
      }
    }
    if(iNoiseTest){
      phz[iOutIDX] = dhzavg / (double)iChans;
      psat[iOutIDX] = dsatavg / (double)iChans;
      if(pdelta) pdelta[iOutIDX] /= (double)iChans;
      if(ptheta) ptheta[iOutIDX] /= (double)iChans;
      if(pbeta) pbeta[iOutIDX] /= (double)iChans;
      if(pgamma) pgamma[iOutIDX] /= (double)iChans;
      if(pripple) pripple[iOutIDX] /= (double)iChans;
    }

    
    int iC1,iC2;
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=0;iC2<=iC1;iC2++){
	if(iC1==iC2){
	  pCorrel[iC1][iC2]=1.0; 
	  continue;
	} else pCorrel[iC1][iC2]=0.; 
	for(i=0;i<iWinSz;i++){
	  pCorrel[iC1][iC2] += pTSN[iC1][i]*pTSN[iC2][i];
	}
	pCorrel[iC1][iC2] /= (double)iWinSz;
	pCorrel[iC2][iC1] = pCorrel[iC1][iC2];
      }
    }
    double dpsum = 0.0; int icc = 0; 
    if(idoabs){
      for(iC1=0;iC1<iChans;iC1++){
	for(iC2=0;iC2<iC1;iC2++){
	  dpsum += fabs(pCorrel[iC1][iC2]);
	  icc++;
	}
      }
    } else {
      for(iC1=0;iC1<iChans;iC1++){
	for(iC2=0;iC2<iC1;iC2++){
	  dpsum += pCorrel[iC1][iC2];
	  icc++;
	}
      }
    }
    dpsum /= (double) icc;
    if(dpsum > 1. || dpsum < -1.){
      printf("out of bounds = %g!\n",dpsum);
    }
    vTCor[iOutIDX] = dpsum;
  }
  dRet = 1.0;
  printf("\n");

CLEANUP:
  if(pTS) FreeListVec(&pTS);
  if(pCorrel) freedouble2D(&pCorrel,iChans);
  if(pTSN) freedouble2D(&pTSN,iChans);
  free(psums);
  free(psums2);
  if(pfft) free(pfft);
  return dRet;
  ENDVERBATIM
}


FUNCTION tscor () {
  VERBATIM

  int i;

  double dRet = 0.0;

  ListVec* pEIG = 0, *pTS = 0;

  double** pTSN = 0; 
  double** pCorrel = 0; 
  double** pV = 0;
  double* pW = 0; 

  pEIG = AllocListVec(*hoc_objgetarg(1)); 
  pTS = AllocListVec(*hoc_objgetarg(2)); 
  
  if(!pEIG || !pTS){
    printf("tscor ERRA: problem initializing 1st 2 args!\n");
    goto CLEANUP;
  }

  int iChans = pTS->isz; 
  if(iChans < 2) {
    printf("tscor ERRB: must have at least 2 EEG channels!\n");
    goto CLEANUP;
  }
  int iWinSz = (int) *getarg(3); 
  int iINC = (int) *getarg(4); 

  int iSamples = pTS->plen[0]; 

  if(verbose) printf("iWinSz=%d iINC=%d iSamples=%d\n",iWinSz,iINC,iSamples);
  for(i=1;i<iChans;i++) {
    if(pTS->plen[i] != iSamples){
      printf("tscor ERRC: time series of unequal size %d %d!\n",iSamples,pTS->plen[i]);
      goto CLEANUP;
    }
  }
  int iOutSz = (iSamples - iWinSz + 1) / iINC;
  if(pEIG->isz < iOutSz){
    printf("tscor ERRD: need List of size at least %d as arg 1, have only %d!\n",iOutSz,pEIG->isz);
    goto CLEANUP;
  }
  if(verbose) printf("iOutSz=%d\n",iOutSz);

  pCorrel = getdouble2D(iChans,iChans); if(verbose) printf("got pCorrel\n");
  pTSN = getdouble2D(iChans,iWinSz); if(verbose) printf("got pTSN\n");
  pV = getdouble2D(iChans,iChans);   if(verbose) printf("got pV\n");
  pW = (double*) malloc(sizeof(double)*iChans); if(verbose) printf("got pW\n");
  if(!pCorrel || !pTSN || !pV || !pW){
    printf("tscor ERRD: out of memory!\n");
    goto CLEANUP;
  }

  double dVal = 0.0;
  int iStartIDX = 0 , ierr = 0 , iTIDX = 0;
  for(iStartIDX=0;iStartIDX+iWinSz<iSamples;iStartIDX+=iINC,iTIDX++){
    if(iStartIDX%100==0) printf("%d\n",iStartIDX);
    int IDX = iStartIDX , iChan = 0 , iEndIDX = iStartIDX + iWinSz;
    int iCurSz = iEndIDX - iStartIDX;
    for(iChan=0;iChan<iChans;iChan++){ 
      double dSum = 0.0, dSum2 = 0.0;
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++){
        dVal = pTS->pv[iChan][IDX];
	dSum += dVal;
	dSum2 += dVal*dVal;
      }
      dSum /= iCurSz; 
      dSum2 /= iCurSz; 
      dSum2 -= dSum*dSum;
      dSum2 = sqrt(dSum2); 
      i=0;
      for(IDX=iStartIDX;IDX<iEndIDX;IDX++,i++){
	pTSN[iChan][i] = (pTS->pv[iChan][IDX] - dSum) / dSum2;
      }
    }
    
    int iC1,iC2;
    for(iC1=0;iC1<iChans;iC1++){
      for(iC2=0;iC2<=iC1;iC2++){
	if(iC1==iC2){
	  pCorrel[iC1][iC2]=1.0; 
	  continue;
	} else pCorrel[iC1][iC2]=0.;
	for(i=0;i<iWinSz;i++){
	  pCorrel[iC1][iC2] += pTSN[iC1][i]*pTSN[iC2][i];
	}
	pCorrel[iC1][iC2] /= (double)iWinSz;
	pCorrel[iC2][iC1] = pCorrel[iC1][iC2];
      }
    }
    
    mysvd(iChans,iChans,pCorrel,pW,pV,&ierr); 
    qsort(pW,iChans,sizeof(double),mycompare); 

    memcpy(pEIG->pv[iTIDX],pW,sizeof(double)*iChans); 
  }
  dRet = 1.0;
  printf("\n");

CLEANUP:
  if(pEIG) FreeListVec(&pEIG);
  if(pTS) FreeListVec(&pTS);
  if(pCorrel) freedouble2D(&pCorrel,iChans);
  if(pTSN) freedouble2D(&pTSN,iChans);
  if(pW) free(pW);
  if(pV) freedouble2D(&pV,iChans);
  return dRet;
  ENDVERBATIM
}

VERBATIM

double getmean(double* p,int n) {
  if(n<1) return 0.0;
  int i = 0;
  double sum = 0.0, *pp=p;
  for(i=0;i<n;i++,pp++) sum += *pp;
  return sum / (double) n;
}

double getstd(double* p,int n,double X) {
  if(n<1) return 0.;
  int i;
  double X2=0., *pp=p;
  for(i=0;i<n;i++,pp++) X2 += (*pp * *pp);
  X2 = X2/(double)n - X*X;
  if(X2>0.) return sqrt(X2);
  return 0.;
}

double dnormv(double* p,int n) {
  if(n<1) return 0.0;
  double X = getmean(p,n);
  double s = getstd(p,n,X);
  double* pp = p;
  int i;
  if(s>0.)
    for(i=0;i<n;i++,pp++) *pp = (*pp - X)/s; 
  else
    for(i=0;i<n;i++,pp++) *pp -= X;
  return 1.0;
}

static double normv(void* v){ 
  double* x;
  int sz = vector_instance_px(v,&x);
  if(sz<1){
    printf("normv ERRA: empty input size!\n");
    return 0.;
  }
  return dnormv(x,sz);
}

double dcopynzidx(double* pin,double* pidx,double* pout,int sz) {
  int szout,i;
  szout=0;
  for(i=0;i<sz;i++) if(pidx[i]) pout[szout++]=pin[i];
  return (double) szout;
}

static double copynzidx (void* v) {
  double *x,*y,*z,ret; int sz;
  sz = vector_instance_px(v,&x);
  if(sz!=vector_arg_px(1,&y)) y=vector_newsize(vector_arg(1),sz);
  if(sz!=vector_arg_px(2,&z)) z=vector_newsize(vector_arg(2),sz);
  ret=dcopynzidx(x,y,z,sz);
  vector_resize(vector_arg(2),(int)ret);
  return ret;
}
ENDVERBATIM


FUNCTION tsmul () {
  VERBATIM
  double* x, *y, *px, *py, sum = 0.;
  int xsz,ysz,i;
  xsz = vector_arg_px(1,&x);
  ysz = vector_arg_px(2,&y);
  px=x; 
  py=y;
  if(xsz!=ysz || xsz<1 || ysz<1){
    printf("tsmul ERRB: input vecs are invalid sizes: %d %d!\n",xsz,ysz);
    return -2.;
  }
  for(i=0;i<xsz;i++,px++,py++) sum += *px * *py;
  return sum / (double) xsz;
  ENDVERBATIM
}

PROCEDURE install () {
  if (INSTALLED==1) {
    printf("already installed tsa.mod")
  } else {
    INSTALLED=1
    VERBATIM
      
    install_vector_method("inrange",inrange);
    install_vector_method("normv",normv);
    install_vector_method("copynzidx",copynzidx);
    ENDVERBATIM
  }
}