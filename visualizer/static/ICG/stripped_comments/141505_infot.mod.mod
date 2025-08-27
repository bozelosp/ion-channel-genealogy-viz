NEURON {
  SUFFIX infot
  GLOBAL installed,beg,end
  GLOBAL verbose
  GLOBAL MINEXP,MAXEXP,MINLOG2,MAXLOG2,count,cutoff,binmin,binmax
}

PARAMETER {
  installed = 0
  verbose = 0.5
  useslice = 0 
  beg = 0 
  end = 0
  MINEXP=-20
  MAXEXP=20
  MINLOG2=0.001
  MAXLOG2=32
  count=0
  cutoff=0.2 
  binmin=0 
  binmax=0 
  KTProb=0 
}

VERBATIM
#include "misc.h"
#include <stdlib.h>
#include <math.h>
static const double* ITsortdata = NULL; 
static double tetrospks2(), pdfpr(), tetrospks3();
static int dbxi[10];

typedef struct ITNode_ {
  int idims;
  int icount;
  int* pvals;
  struct ITNode_* pLeft;
  struct ITNode_* pRight;
} ITNode;

ITNode* allocITNode(int idims) {
  ITNode* p;
  p = calloc(1,sizeof(ITNode));
  if(!p) { printf("allocITNode: out of mem!\n"); hxe(); return 0x0; }
  p->pvals = calloc(idims,sizeof(int));
  return p;
}

void freeITNode(ITNode** pp) {
  ITNode* p;
  p = pp[0];
  free(p->pvals);
  if(p->pLeft != 0x0) freeITNode(&p->pLeft);
  if(p->pRight!= 0x0) freeITNode(&p->pRight);
  free(p);
  *pp = 0x0;
}

int compITNode(ITNode* pL,ITNode* pR) {
  int idims,i;
  idims = pL->idims;
  for(i=0;i<idims;i++) {
    if(pL->pvals[i] < pR->pvals[i]) return -1; 
    if(pL->pvals[i] > pR->pvals[i]) return 1;  
  }
  return 0; 
}



ITNode* findITNode( ITNode* pRoot, ITNode* pS ) {
  int ic;
  if(pRoot == 0x0 || pS == 0x0) return 0x0;
  ic = compITNode(pS,pRoot);
  if(ic==0) {
    return pRoot;
  } else if(ic == -1) {
    return findITNode(pRoot->pLeft,pS);
  } else return findITNode(pRoot->pRight,pS);
}





int insertITNode( ITNode* pRoot, ITNode* pNew ) {
  int ic;
  ic = compITNode(pNew,pRoot);
  if(ic == 0) {
    pRoot->icount++;
    return 0;
  }
  if(ic == -1) {
    if(pRoot->pLeft == 0x0) {
      pRoot->pLeft = pNew;
      return 1;
    } 
    return insertITNode(pRoot->pLeft,pNew);
  } else {
    if(pRoot->pRight == 0x0) {
      pRoot->pRight = pNew;
      return 1;
    }
    return insertITNode(pRoot->pRight,pNew);
  }
}

ENDVERBATIM

FUNCTION TestITree () {
  VERBATIM
  int i,j;
  ITNode *pRoot,pNew;
  return 0.0;
  ENDVERBATIM
}

VERBATIM


double getavgd(double* p,int n) {
  int i; double sum,*pp;
  if(n<1) return 0.0;
  sum = 0.0; pp=p;
  for(i=0;i<n;i++,pp++) sum += *pp;
  return sum / (double) n;
}

double getstdevd(double* p,int n) {
  int i;
  double X,X2,*pp;
  if(n<1) return 0.;
  X=X2=0.0; pp=p;
  for(i=0;i<n;i++,pp++) {
    X2 += (*pp * *pp);
    X += *pp;
  }
  X /= (double)n;
  X2 = X2/(double)n - X*X;
  if(X2>0.) return sqrt(X2);
  return 0.;
}

double getstdevi(int* p,int n) {
  double sd,X,X2;
  int i;
  X = X2 = 0.0;
  for(i=0;i<n;i++) {
    X += p[i];
    X2 += p[i]*p[i];
  }
  X /= (double) n;
  sd = X2/n - X*X;
  if(sd>0.) return sqrt(sd);
  return 0.0;
}

double kprob1Dd (double* X,int sz,double h,double val) {
  double prob,dif,sig2,h2;
  int i;
  h2 = 2.0*h*h;
  prob = 0.0;
  for(i=0;i<sz;i++) {
    dif = val-X[i];
    dif *= dif;
    prob += exp( -dif / h2 );
  }
  prob *= (1.0 / ( (double) sz * h * SQRT2PI ) );
  return prob;
}

double kprob2D (int* X,int* Y,int sz,double sigma,double a,int x,int y) {
  double prob,dx,dy,sig2; 
  int i;
  prob = 0.0;
  sig2 = sigma * sigma;
  for(i=0;i<sz;i++) {
    dx = x - X[i];
    dy = y - Y[i];
    dx *= dx;
    dy *= dy;
    prob += exp( (-dx-dy) / sig2 );
  }
  return prob;
}

double log2d ( double d ) {
  return log(d) / LG2;
}





double kprobstepi (int** XX,int iDims,int iStartIDX,int iLen,double h,int* x) {
  double prob,dif,val,cnt;
  int i,j;
  prob = cnt = 0.0;
  for(i=iStartIDX;i<iLen;i++) { 
    dif = 0.0;
    for(j=0;j<iDims;j++) { 
      val = x[j] - XX[i][j];
      dif += val*val;
    }
    dif = h - sqrt(dif);
    if(dif > 0) prob += 1.0; 
    cnt += 1.0;
  }  
  prob = prob / cnt; 
  return prob;
}





double kprobgaussi (int** XX,int iDims,int iStartIDX,int iLen,double h,int* x) {
  double prob,dif,val,h2,hh,tp,tpn;  
  int i,j;
  tp = 2.0*3.14159265;
  tpn = pow(tp,((double)iDims)/2.0);
  h2 = h; 
  hh = h;



  prob = 0.0;
  for(i=iStartIDX;i<iLen;i++) {
    dif = 0.0;
    for(j=0;j<iDims;j++) {
      val = x[j] - XX[i][j];
      dif += val*val;
    }
    prob += exp( -dif / h2 );
  }

  prob = prob / ( tpn * iLen );
  return prob;
}


double getbandwidthd (int d,int N,double sd) {

  
  if(d > 1) {
    return sd*pow((4.0/(double)d+2.0),1.0/((double)d+4.0))*pow((double)N,-1.0/(d+4.0));
  } else {
    return sd * 1.06 * pow((double)N,-0.2);
  }

  
  

  
  

  return 1.;
}


int* getnormd (double* x,int sz,int nbins) {
  int* p;
  double max,min,rng,dnbins;
  int i;
  max=min=x[0];
  for(i=1;i<sz;i++) {
    if(x[i]>max) max=x[i];
    if(x[i]<min) min=x[i];
  }
  p = (int*) malloc(sizeof(int)*sz);
  rng = max - min; dnbins = nbins;
  for(i=0;i<sz;i++) {
    p[i] = (int) ((dnbins * ((double)x[i]-min)/rng) + 0.5);
    if(p[i] < 0) p[i] = 0; else if(p[i] >= nbins) p[i] = nbins - 1;
  }
  return p;
}



int downsampavgd (double* pin,double* pout,int szin,int szout,int winsz) {
  double dsz,dsum;
  dsz= (double)szin/(double)winsz;  dsz = ceil(dsz);
  if(szout<dsz) return 0;
  int i,j,k,startx,endx,cnt;
  for(i=0,k=0;i<szin;i+=winsz) {
    startx=i; endx=i+winsz-1; if(endx>=szin) endx=szin-1;
    cnt=endx-startx+1; dsum=0.;
    for(j=startx;j<=endx;j++) dsum += pin[j];
    if(cnt>0) dsum /= (double) cnt;
    pout[k++] = dsum;
  }
  return 1;
}


void chgencd (double* pin,double* pout,int sz) {
  double df; int i;
  if(sz<1) return;
  pout[0]=1;
  for(i=1;i<sz;i++) {
    df=pin[i]-pin[i-1];
    if(df>0) {
      pout[i]=2;
    } else if(df<0) {
      pout[i]=0;
    } else pout[i]=1;
  }
}


int oddometerinc (int* p,int idims,int imaxval) {
  int i;
  i = 0;
  for(;i<idims;i++) {
    if(p[i] + 1 > imaxval) {
      p[i] = 0;
      continue;
    }
    p[i] += 1;
    return 1;
  }
  return 0;
}

int* doublep2intp (double* d, int sz) {
  int i;
  int* pout;
  pout = (int*) malloc(sizeof(int)*sz);
  if(!pout) { printf("doublep2intp: out of mem!\n"); hxe(); }
  for(i=0;i<sz;i++) pout[i] = (int) d[i]; 
  return pout;
}

void pri2d (int** p, int r,int c) {
  int i,j;
  for(i=0;i<r;i++) {
    for(j=0;j<c;j++) {
      printf("%d ",p[i][j]);
    }
    printf("\n");
  }
}




int issubint (int* str,int sz1,int* pat,int sz2) {
  int i,j;  
  for(i=0;i+sz2<=sz1;i++) {
    for(j=0;j<sz2 && str[i+j]==pat[j];j++);
    if(j==sz2) return 1;
  }
  return 0;
}



double lz76complexityd (double* x,int iLen) {
  int *px,**pdict,i,j,npat,patlen,fnd;
  if(iLen<1) return 0.0;
  px=doublep2intp(x,iLen);
  pdict=getint2D(iLen,2);
  npat=i=1; pdict[0][0]=0; pdict[0][1]=1;
  while(i<iLen) {    
    fnd=1;
    for(patlen=1;i+patlen<iLen;patlen++) { 
      if(verbose>2) printf("i%d patlen%d\n",i,patlen);
      if(!(fnd=issubint(px,i,&px[i],patlen))){
        pdict[npat][0]=i; pdict[npat][1]=patlen; npat++; i+=patlen; if(verbose>2) printf("not found\n"); break;
      } 
    }
    if(fnd) {pdict[npat][0]=i; pdict[npat][1]=iLen-i; npat++; break;} 
  }
  if(verbose>1) { 
    j=0; printf("found %d patterns:\n\t",npat);
    for(i=0;i<iLen;i++) {
      if(i==pdict[j][0]+pdict[j][1]){ printf("|"); j++; }
      printf("%d",px[i]);
    }  printf("\n");
    printf("pat off len:\n");
    for(i=0;i<npat;i++) printf("\t%d %d\n",pdict[i][0],pdict[i][1]);
  }
  free(px); freeint2D(&pdict,iLen);
  return log2d((double)iLen) * npat / (double)iLen;
}


static double lz76c (void* vv) {
  double* x; int n;
  n = vector_instance_px(vv,&x);
  if(n<1) { printf("lz76c ERRA: 0 size vector\n"); return 0.0; }
  return lz76complexityd(x,n);
}








double tentropd (double* x,double* y,int iLen,int nbins,int xpast,int ypast,int nshuf,double* sig,double hval,int kern,double* nTE) {
  int *px,*py,jointd,i,j,k,**pjointxy,minidx,**pjointx,**ppastx,**ppastxy,*poddometer,sh,*ptmp;
  int cntjxy,cntjx,cntpx,cntpxy;
  double teout,te,teavg,testd,probjxy,probjx,probpx,probpxy;
  double hjxy,hjx,hpx,hpxy,sd,dsum,dtmp,N,ent,norm;
  px=py=poddometer=0x0; pjointxy=pjointx=ppastx=ppastxy=0x0;
  dsum=0.0;
  *sig=-1e6;
  *nTE = 0.0;
  sh=cntjxy=cntjx=cntpx=cntpx=0;   teavg=testd=te=teout=0.0;

  px = nbins>0 ? getnormd(x,iLen,nbins) : doublep2intp(x,iLen); 
  py = nbins>0 ? getnormd(y,iLen,nbins) : doublep2intp(y,iLen); 

  jointd = xpast + ypast + 1;

  if(verbose>1) printf("iLen%d, jointd%d, xpast%d, ypast%d\n", iLen,jointd,xpast,ypast);

  
  pjointxy=getint2D(iLen,jointd+(kern==0?1:0)); if(!pjointxy) {printf("tentropd: out of mem!\n"); hxe();}
  pjointx=getint2D(iLen,xpast+1+(kern==0?1:0)); if(!pjointx)  {printf("tentropd: out of mem!\n"); hxe();}

  ptmp=kern?0x0:(int*)malloc(sizeof(int)*jointd);

  if(verbose>1) printf("jointd=%d, xpast+1=%d, iLen=%d, nbins=%d\n",jointd,xpast+1,iLen,nbins);

  if(xpast>0) ppastx=getint2D(iLen,xpast+(kern==0?1:0)); else {ppastx=0; printf("tentropd WARN: ppastx=0!\n");}
  if(xpast>0 || ypast>0) {
    ppastxy=getint2D(iLen,xpast+ypast+(kern==0?1:0));
  } else {ppastxy=0;  printf("tentropd WARN: ppastxy=0!\n");}

  minidx = xpast > ypast ? xpast : ypast;

  do {

    if(sh > 0) ishuffle(py,iLen); 

    if(kern) { 
      for(i=minidx;i<iLen-1;i++) { 
        pjointxy[i][0] = pjointx[i][0] = px[i+1]; 
        k = 1;
        for(j=0;j<xpast;j++,k++) {
          pjointxy[i][k] = pjointx[i][k] = ppastx[i][k-1] = ppastxy[i][k-1] = px[i-j]; 
        }
        for(j=0;j<ypast;j++,k++) {
          pjointxy[i][k] = ppastxy[i][k-1] = py[i-j]; 
        }
      }
      if(sh==0){ 
        if(hval<=0.0){   
          sd = (getstdevi(px,iLen) + getstdevi(py,iLen)) / 2.0;
          if(verbose>1) printf("sd=%g\n",sd);
          hjxy = getbandwidthd(jointd,iLen-1-minidx+1,sd);
          hjx = getbandwidthd(1+xpast,iLen-1-minidx+1,sd);
          hpx = getbandwidthd(xpast,iLen-1-minidx+1,sd);
          hpxy = getbandwidthd(xpast+ypast,iLen-1-minidx+1,sd);
        } else hjxy=hjx=hpx=hpxy=hval;
        if(verbose>1) printf("hjxy=%g, hjx=%g, hpx=%g, hpxy=%g\n",hjxy,hjx,hpx,hpxy);
      }    
      
      te=dsum=0.0;
      for(j=minidx;j<iLen-1;j++) {
        poddometer = pjointxy[j];
        probjxy = kprobstepi(pjointxy,jointd,minidx,iLen-1,hjxy,poddometer);
        dsum += probjxy;
        probjx = kprobstepi(pjointx,xpast+1,minidx,iLen-1,hjx,poddometer);
        probpx = xpast>0?kprobstepi(ppastx,xpast,minidx,iLen-1,hpx,&poddometer[1]):1.0;
        probpxy = kprobstepi(ppastxy,xpast+ypast,minidx,iLen-1,hpxy,&poddometer[1]);    
        if(verbose>2) printf("probjxy=%g,probjx=%g,probpx=%g,probpxy=%g\n",probjxy,probjx,probpx,probpxy);
        if(probpx>1e-9 && probjxy>1e-9 && probjx>1e-9 && probpxy>1e-9) {
          dtmp = (probpx*probjxy)/(probjx*probpxy);
          if(usetable && dtmp>=MINLOG2 && dtmp<=MAXLOG2) 
           te += probjxy * _n_LOG2(dtmp);
          else
           te += probjxy * log2d(dtmp);
        }
      }
      N = ( iLen - 1 ) - minidx + 1;
      if(N > 0) dsum = dsum / (double) N;
      if(verbose>1) printf("dsum=%g\n",dsum);
    } else { 

      N = ( iLen - 1 ) - minidx + 1;

      
      if(sh==0) cntjxy=cntjx=cntpx=cntpxy=0; else cntjxy=cntpxy=0;

      for(i=minidx;i<iLen-1;i++) { 
        ptmp[0]=px[i+1]; 
        k=1;
        for(j=0;j<xpast;j++,k++) ptmp[k] = px[i-j]; 
        for(j=0;j<ypast;j++,k++) ptmp[k] = py[i-j]; 
        for(j=0;j<cntjxy;j++) { 
          for(k=0;k<jointd;k++) if(ptmp[k]!=pjointxy[j][k])break;
          if(k==jointd)break;
        }
        if(j==cntjxy) { 
          for(k=0;k<jointd;k++) pjointxy[cntjxy][k]=ptmp[k];
          pjointxy[cntjxy][k] = 1; 
          cntjxy++; 
        } else {
          pjointxy[j][jointd] += 1; 
        }
      }
      if(verbose>1) {printf("cntjxy=%d\n",cntjxy);
        for(i=0;i<cntjxy;i++) {
          printf("pjointxy%d: [");
          for(j=0;j<jointd;j++) printf("%d ",pjointxy[i][j]);
          printf("] , cnt=%d, p=%g\n",pjointxy[i][jointd],pjointxy[i][jointd]/(double)N);
        }
      }
      if(sh==0) { 
        for(i=minidx;i<iLen-1;i++) { 
          ptmp[0]=px[i+1];
          k=1;
          for(j=0;j<xpast;j++,k++) ptmp[k] = px[i-j]; 
          for(j=0;j<cntjx;j++) { 
            for(k=0;k<xpast+1;k++) if(ptmp[k]!=pjointx[j][k])break;
            if(k==xpast+1)break;
          }
          if(j==cntjx) { 
            for(k=0;k<xpast+1;k++) pjointx[cntjx][k]=ptmp[k];
            pjointx[cntjx][k] = 1; 
            cntjx++; 
          } else {
            pjointx[j][xpast+1] += 1; 
          }        
        }
        for(i=minidx;i<iLen-1;i++) { 
          k=1;
          for(j=0;j<xpast;j++,k++) ptmp[k-1] = px[i-j]; 
          for(j=0;j<cntpx;j++) { 
            for(k=0;k<xpast;k++) if(ptmp[k]!=ppastx[j][k])break;
            if(k==xpast)break;
          }
          if(j==cntpx) { 
            for(k=0;k<xpast;k++) ppastx[cntpx][k]=ptmp[k];
            ppastx[cntpx][k] = 1; 
            cntpx++; 
          } else {
            ppastx[j][xpast] += 1; 
          }                
        }
      }
      for(i=minidx;i<iLen-1;i++) { 
        k = 1;
        for(j=0;j<xpast;j++,k++) ptmp[k-1] = px[i-j]; 
        for(j=0;j<ypast;j++,k++) ptmp[k-1] = py[i-j]; 
        for(j=0;j<cntpxy;j++) { 
          for(k=0;k<xpast+ypast;k++) if(ptmp[k]!=ppastxy[j][k])break;
          if(k==xpast+ypast)break;
        }
        if(j==cntpxy) { 
          for(k=0;k<xpast+ypast;k++) ppastxy[cntpxy][k]=ptmp[k];
          ppastxy[cntpxy][k] = 1; 
          cntpxy++; 
        } else {
          ppastxy[j][xpast+ypast] += 1; 
        }                        
      }
      te=dsum=0.0;
      if(sh==0) { 
        norm = 0.0;
        for(i=0;i<cntjx;i++) {
          probjx = pjointx[i][xpast+1]/(double)N;
          for(j=0;j<cntpx;j++) { 
            for(k=0;k<xpast;k++) { if(ppastx[j][k]!=pjointx[i][k+1]) break; }
            if(k==xpast) {
              probpx = ppastx[j][xpast]/(double)N;
              break;
            }
          }
          if(probjx > 1e-9 && probpx > 1e-9) {
            dtmp = probjx / probpx;
            if(usetable && dtmp>=MINLOG2 && dtmp<=MAXLOG2) 
              norm -= probjx * _n_LOG2(dtmp);
            else
              norm -= probjx * log2d(dtmp);            
          }
        }      
      }
      ent = 0.0;
      for(i=0;i<cntjxy;i++) { 

        probjxy = pjointxy[i][jointd]/(double)N; 
        dsum += probjxy;
        
        for(j=0;j<cntpxy;j++) { 
          for(k=0;k<xpast+ypast;k++) { if(ppastxy[j][k]!=pjointxy[i][k+1]) break; }
          if(k==xpast+ypast) {
            probpxy = ppastxy[j][xpast+ypast]/(double)N;
            break;
          }
        }
        if(verbose>1) printf("probjxy=%g,probjx=%g,probpx=%g,probpxy=%g\n",probjxy,probjx,probpx,probpxy);
        if(probjxy>1e-9 && probpxy>1e-9) {
          dtmp = probjxy / probpxy;
          if(usetable && dtmp>=MINLOG2 && dtmp<=MAXLOG2) 
            ent -= probjxy * _n_LOG2(dtmp);
          else
            ent -= probjxy * log2d(dtmp);
        }
      }
      te = norm - ent;
      if(te<0) te=0.0; 
      if(verbose>1) printf("dsum=%g\n",dsum);
    }
    if(sh==0) teout=te; else {
      teavg += te;
      testd += te*te;
    }
  } while(sh++ < nshuf);

  if(nshuf>0) {
    teavg /= (double) nshuf;
    testd = testd/(double)nshuf - teavg*teavg;
    if(testd > 1e-9) *sig = (teout - teavg) / testd;
    if(verbose>0) printf("teout=%g,teavg=%g,testd=%g,sig=%g\n",teout,teavg,testd,*sig);
    *nTE = norm > 1e-9 ? (teout - teavg) / norm : (teout - teavg);
  } else *nTE=norm>1e-9?teout/norm:teout;
  if(*nTE < 0) *nTE = 0.; 
  freeint2D(&pjointxy,iLen); freeint2D(&pjointx,iLen);   
  if(ppastx) freeint2D(&ppastx,iLen); if(ppastxy) freeint2D(&ppastxy,iLen);  
  free(px); free(py); if(ptmp) free(ptmp);
  return teout;
}








static double tentrop (void* vv) {
  double *x,*y,hval,sig,*pout,te,nTE;
  int szx,szy,szo,nbins,xpast,ypast,nshuf,kern;
  szx = vector_instance_px(vv,&x);
  if((szy=vector_arg_px(1,&y))!=szx) {
    printf("tentrop ERRA: x,y must have same size (%d,%d)\n",szy,szy);
    return -1.0;
  } 
  
  nbins = *getarg(2);
  xpast = ifarg(3)?(int)*getarg(3):1;
  ypast = ifarg(4)?(int)*getarg(4):2;
  nshuf = ifarg(5)?(int)*getarg(5):0;
  if(ifarg(6)) szo=vector_arg_px(6,&pout); else szo=0; 
  hval = ifarg(7)?*getarg(7):-1.0;
  kern = ifarg(8)?(int)*getarg(8):0;
  te = tentropd(y,x,szx,nbins,ypast,xpast,nshuf,&sig,hval,kern,&nTE);
  if(szo>0) pout[0]=te; if(szo>1) pout[1]=sig;
  return nTE;
}


void myprvec (double* x, int sz,char* msg) {
  int i;
  if(msg) printf("%s",msg);
  for(i=0;i<sz;i++) printf("%g ",x[i]);
  printf("\n");
}











static double ntedir (void* vv) {
  double *x,*y,hval,sigxy,sigyx,*pout,texy,teyx,texyout,teyxout,prefdir;
  double ntexy,nteyx,ntexyout,nteyxout;
  int szx,szy,szo,nbins,xpast,ypast,nshuf,i,kern;
  szx = vector_instance_px(vv,&x);
  if((szy=vector_arg_px(1,&y))!=szx) {
    printf("tecause ERRA: x,y must have same size (%d,%d)\n",szx,szy);
    return -1.0;
  } 
  
  
  nbins = *getarg(2);
  xpast = ifarg(3)?(int)*getarg(3):1;
  ypast = ifarg(4)?(int)*getarg(4):2;
  nshuf = ifarg(5)?(int)*getarg(5):0;
  if(ifarg(6)) szo=vector_arg_px(6,&pout); else szo=0; 
  hval = ifarg(7)?*getarg(7):-1.0;
  kern = ifarg(8)?(int)*getarg(8):0;

  texyout = tentropd(y,x,szx,nbins,ypast,xpast,nshuf,&sigxy,hval,kern,&ntexyout); 
  teyxout = tentropd(x,y,szy,nbins,ypast,xpast,nshuf,&sigyx,hval,kern,&nteyxout); 

  if(verbose) printf("texyout=%g,teyxout=%g,ntexyout=%g,nteyxout=%g\n",texyout,teyxout,ntexyout,nteyxout);

  prefdir = ntexyout + nteyxout > 1e-9 ? (ntexyout - nteyxout)/(ntexyout+nteyxout) : 0;

  if(szo>0) pout[0]=texyout; if(szo>1) pout[1]=sigxy; if(szo>2) pout[2]=ntexyout;
  if(szo>3) pout[3]=teyxout; if(szo>4) pout[4]=sigyx; if(szo>5) pout[5]=nteyxout;
  if(szo>6) pout[6]=prefdir;

  return prefdir;
}








static double tentropspks (void* vv) {
  double *X1,*X2,*XO,*X3; int szX1,szX2,szXO,shuf,szX3;
  szX1 = vector_instance_px(vv,&X1);
  if((szX2=vector_arg_px(1,&X2))!=szX1) {
    printf("tentropspks ERRA: X1,X2 must have same size (%d,%d)\n",szX1,szX2); return -1.0; }
  szXO=ifarg(2)?vector_arg_px(2,&XO):0;
  shuf=ifarg(3)?((int)*getarg(3)):0;
  szX3=ifarg(4)?vector_arg_px(4,&X3):0;
  if(szX3) tetrospks3(X1,X2,X3,XO,szX1,szXO,shuf);
  else return tetrospks2(X1,X2,XO,szX1,szXO,shuf);
}


static double ntel2 (void* vv) {
  double *pfrom,*pto,*pTE,*pNTE,ret,te;
  int i,j,sz,tsz,nshuf;
  ListVec *pLFrom,*pLTo;
  ret = 0.0; 
  sz = vector_instance_px(vv, &pfrom);
  pTE=0x0;
  if(sz!=vector_arg_px(1, &pto) ||  sz!=vector_arg_px(4,&pNTE) || (ifarg(5) && sz!=vector_arg_px(5,&pTE)) ) {
    printf("ntel2 ERRA: arg 1,4,5 must have size of %d\n",sz); hxe();}
  nshuf = ifarg(6) ? (int)*getarg(6):30;
  pLFrom = AllocListVec(*hoc_objgetarg(2)); 
  pLTo = AllocListVec(*hoc_objgetarg(3)); 
  if(!pLFrom->isz || !pLTo->isz) { printf("ntel2 ERRB: empty list vectors\n"); goto NTEL2DOFREE; }
  tsz = pLFrom->plen[0]; 
  if(pTE) {
    for(i=0;i<sz;i++,pfrom++,pto++,pNTE++,pTE++) {
      *pNTE=tetrospks2(pLFrom->pv[(int)*pfrom],pLTo->pv[(int)*pto],&te,tsz,1,nshuf);
      *pTE=te;
    }
  } else {
    for(i=0,j=0;i<sz;i++,pfrom++,pto++,pNTE++) {
      *pNTE=tetrospks2(pLFrom->pv[(int)*pfrom],pLTo->pv[(int)*pto],0x0,tsz,0,nshuf);
    }
  }
  ret=1.0;
NTEL2DOFREE:
  FreeListVec(&pLFrom); FreeListVec(&pLTo);
  return ret;
}


static double nte (void* vv) {
  int i, j, k, ii, jj, oi, nx, ny, omax, sz, ci, bg, flag, nshuf, szX, *x1i, *x2i, x1z, x2z;
  ListVec *pLi; Object *obi,*ob;  double *x, *y, *out, o1, o2, cnt, *vvo[3]; void *vvl;
  nx = vector_instance_px(vv, &x); 
  ny = vector_arg_px(1, &y);       
  pLi = AllocListVec(obi=*hoc_objgetarg(2)); 
  ob =   *hoc_objgetarg(3);
  if (ISVEC(ob)) {omax  = vector_arg_px(3, &out);} else omax=-1;
  if (ifarg(4)) nshuf=(int)*getarg(4); else nshuf=0;
  if (nshuf<=0) nshuf=0;
  if (nshuf>200) {printf("WARN: reducing # of shuffles from %d to 200\n",nshuf); nshuf=200;}
  if (omax==-1) { flag=1; 
    i = ivoc_list_count(ob);
    if (i<3) {printf("nte() ERRA out list should be at least 3 (%d)\n",i); hxe();}
    for (k=0;k<3;k++) { 
      j=list_vector_px3(ob, k, &vvo[k], &vvl);
      if (k==0) omax=j; else if (omax!=j||j<10) {
        printf("nte() ERRC: too small %d %d\n",j,omax); hxe(); }
    }
  } else if (2*nx*ny==omax) { flag=4; 
  } else if (nx==ny && 2*nx==omax) { flag=2; 
  } else if (nx==ny && nx==omax) {   flag=3; 
  } else {printf("te_infot ERRA vector sizes are %d %d %d\n",nx,ny,omax);hxe();}
  ci=pLi->isz;
  sz=pLi->plen[0];  
  if (flag==1 || flag==4) { 
    x1i=(int*)calloc(nx,sizeof(int)); x2i=(int*)calloc(ny,sizeof(int)); x1z=x2z=0; 
    bg=(int)beg;
    if (binmin>0 && flag==1) {    
      szX=((end>0.0 && end<=sz)?(int)end:sz) - beg;
      for (i=0;i<nx;i++)   { ii=(int)x[i];
        if (ii<0 || ii>=ci || pLi->plen[ii]!=sz) {
          printf("nte ERRC:%d imax:%d isz:%d sz0:%d\n",ii,ci,pLi->plen[ii],sz); hxe(); }
        for (j=0,cnt=0.;j<szX;j++) if (pLi->pv[ii][j+bg]>0) cnt++;
        if (cnt>(double)szX*binmin && (!binmax || cnt<(double)szX*binmax)) { 
          x1i[x1z]=ii; x1z++; }
      }
      for (i=0;i<ny;i++)   { ii=(int)y[i];
        if (ii<0 || ii>=ci || pLi->plen[ii]!=sz) {
          printf("nte ERRCa i:%d imax:%d isz:%d sz0:%d\n",ii,ci,pLi->plen[ii],sz); hxe(); }
        for (j=0,cnt=0;j<szX;j++) if (pLi->pv[ii][j+bg]>0) cnt++;
        if (cnt>szX*binmin && (!binmax || cnt<szX*binmax)) { x2i[x2z]=ii; x2z++; }
      }
    } else { 
      x1z=nx; x2z=ny;
      for (i=0;i<nx;i++) x1i[i]=(int)x[i]; 
      for (i=0;i<ny;i++) x2i[i]=(int)y[i]; 
    }
    if (flag==4) for (oi=0;oi<omax;oi++) out[oi]=-1;
    for (i=0,oi=0;i<x1z;i++)   { ii=x1i[i];
      for (j=0;j<x2z;j++) { jj=x2i[j]; 
        o1=tetrospks2(pLi->pv[ii],pLi->pv[jj],0x0,sz,0,nshuf);
        if (o1<=-10) { 
          if (verbose>3) printf("tentropspks IGNORING %d:%d %d:%d (%g)\n",i,ii,j,jj,o1);
          continue;
        }
        o2=tetrospks2(pLi->pv[jj],pLi->pv[ii],0x0,sz,0,nshuf);
        if (flag==1) {
          if (oi>=omax-1) { omax*=2; for (k=0;k<3;k++) vvo[k]=list_vector_resize(ob, k, omax); }
          vvo[0][oi]=ii;vvo[1][oi]=jj;vvo[2][oi]=o1;oi++;
          vvo[0][oi]=jj;vvo[1][oi]=ii;vvo[2][oi]=o2;oi++;
        } else { oi=i*ny+j;
          out[oi]=o1;      
          out[oi+nx*ny]=o2;
        }
      }
    }
    if (flag==1) for (k=0;k<3;k++) vvo[k]=list_vector_resize(ob, k, oi);
    return (double)oi;
  } else {
    for (oi=0;oi<omax;oi++) out[oi]=-1;
    for (i=0;i<nx;i++) {
      ii=(int)x[i]; jj=(int)y[i];
      if (ii==jj) { out[i]=0; 
        if (flag==2) out[i+nx]=0;
        continue;
      } 
      if (ii<0 || jj<0 || ii>=ci || jj>=ci || pLi->plen[ii]!=sz || pLi->plen[jj]!=sz) {
        printf("te_infot ERRC bad index or sz i:%d j:%d imax:%d isz:%d jsz:%d sz0:%d\n",\
               ii,jj,ci,pLi->plen[ii],pLi->plen[jj],sz); hxe(); }
      out[i]=tetrospks2(pLi->pv[ii],pLi->pv[jj],0x0,sz,0,nshuf);
      if (flag==2) out[i+nx]= tetrospks2(pLi->pv[jj],pLi->pv[ii],0x0,sz,0,nshuf);
    }
    return out[0];
  }
}




double entropxfgxpd (double* pXP, double* pXFXP,int minv,int maxv,int szp) {
  static double tmp[4];
  int k,l;
  
  for(tmp[0]=0.,k=minv;k<=maxv;k++) for(l=minv;l<=maxv;l++) {
    tmp[1]=pXP[l];  tmp[2]=pXFXP[k*szp+l]; 
    if(tmp[1]>0. && tmp[2]>0.) { tmp[3] = tmp[2]/tmp[1]; 
      if (usetable && tmp[3]>=MINLOG2 && tmp[3]<=MAXLOG2) {
        tmp[0] -=tmp[2]*_n_LOG2(tmp[3]);
      } else {  tmp[0] -=tmp[2]*  log2d(tmp[3]); if (usetable&&verbose>0.4) { 
          printf("WARNA:%g outside of [%g,%g] TABLE\n",tmp[4],MINLOG2,MAXLOG2); }}}}
  return tmp[0];        
}



static double entropxfgxp (void* vv) {
  double *x,*pXP,*pXFXP,dret;
  int sz,minv,maxv,cnt,i,j,szp,*X;
  sz = vector_instance_px(vv,&x);
  cnt=0;
  X=scrset(sz);
  minv=1e9; maxv=-1e9;
  for (i=0;i<sz;i++) { 
    X[i]=(int)x[i]; 
    if (X[i]>0) cnt++;   
    if (X[i]>maxv) maxv=X[i];  if (X[i]<minv) minv=X[i]; 
  }  
  if (minv<0) {
    printf("entropxfgxp ERRA: minimum value must be >= 0:%d\n",minv);hxe();}  
  szp = maxv + 1;
  pXFXP = (double*) calloc(szp*szp,sizeof(double));
  pXP = (double*) calloc(szp,sizeof(double));
  for(i=1;i<sz;i++) pXP[ X[i-1] ]++; 
  for(i=minv;i<=maxv;i++) pXP[i] /= (sz-1);
  if (verbose>2) pdfpr(pXP,szp,1,"pXP");
  for(i=0;i<sz-1;i++) pXFXP[ X[i+1]*szp + X[i] ]++; 
  for(i=minv;i<=maxv;i++) for(j=minv;j<=maxv;j++) pXFXP[i*szp+j]/=(sz-1);
  if (verbose>3) pdfpr(pXFXP,szp,2,"pXFXP");
  dret = entropxfgxpd(pXP,pXFXP,minv,maxv,szp);
  free(pXP); free(pXFXP);
  return dret;
}



double entropx2fgx2px1pd (double* pX2FX2PX1P, double* pX2PX1P, int minv1, int maxv1, int minv2, int maxv2, int szp) {
  static double tmp[4];
  double ent = 0.0;
  int l,k,m;
  for(l=minv2;l<=maxv2;l++) for(k=minv2;k<=maxv2;k++) for(m=minv1;m<=maxv1;m++) {  
    tmp[0]=pX2FX2PX1P[k*szp*szp+l*szp+m]; tmp[1]=pX2PX1P[l*szp+m];
    if (tmp[0]>1e-9 && tmp[1]>1e-9) {
      tmp[2] = tmp[0] / tmp[1];
      if (usetable && tmp[2]>=MINLOG2 && tmp[2]<=MAXLOG2) {
        ent -= tmp[0]*_n_LOG2(tmp[2]);
      } else {    ent -= tmp[0] * log2d(tmp[2]);
        if (usetable&&verbose>0.4) {
          printf("WARNB:%g outside of [%g,%g] TABLE (",tmp[2],MINLOG2,MAXLOG2) ; 
          printf("%g, %g, %g)\n",tmp[0],tmp[1],tmp[2]); }
      }
      if(verbose>2){printf("tmp0=%g, tmp1=%g, tmp2=%g\n",tmp[0],tmp[1],tmp[2]);
        printf("l2d:%g\n",log2d(tmp[2])); printf("ent:%g\n",ent); }
    }
  }
  return ent;
}



double entropx3fgx1px2px3pd (double* pX3FX1PX2PX3P, double* pX1PX2PX3P, 
                            int minv1, int maxv1, int minv2, int maxv2, int minv3, int maxv3, int szp) {
  static double tmp[4];
  double ent = 0.0;
  int l,k,m,n;
  for(l=minv3;l<=maxv3;l++) for(k=minv1;k<=maxv1;k++) for(m=minv2;m<=maxv2;m++) for(n=minv3;n<=maxv3;n++) {  
    tmp[0]=pX3FX1PX2PX3P[l*szp*szp*szp+k*szp*szp+m*szp+n]; tmp[1]=pX1PX2PX3P[k*szp*szp+m*szp+n];
    if (tmp[0]>1e-9 && tmp[1]>1e-9) {
      tmp[2] = tmp[0] / tmp[1];
      if (usetable && tmp[2]>=MINLOG2 && tmp[2]<=MAXLOG2) {
        ent -= tmp[0]*_n_LOG2(tmp[2]);
      } else {    ent -= tmp[0] * log2d(tmp[2]);
        if (usetable&&verbose>0.4) {
          printf("WARNB:%g outside of [%g,%g] TABLE (",tmp[2],MINLOG2,MAXLOG2) ; 
          printf("%g, %g, %g)\n",tmp[0],tmp[1],tmp[2]); }
      }
      if(verbose>2){printf("tmp0=%g, tmp1=%g, tmp2=%g\n",tmp[0],tmp[1],tmp[2]);
        printf("l2d:%g\n",log2d(tmp[2])); printf("ent:%g\n",ent); }
    }
  }
  return ent;
}


double pdfnz (double* pdf, int szp, int dim) {
  double x,ds,cnt; int i,j,k,l,m,n;
  cnt = 0.0;
  if(dim==1) {
    for(m=0;m<szp;m++) if(pdf[m]>0.0) cnt+=1.0;
  } else if(dim==2) {
    for(l=0;l<szp;l++) for(m=0;m<szp;m++) if(pdf[l*szp+m]>0.0) cnt+=1.0;
  } else if(dim==3) {
    for(k=0;k<szp;k++) for(l=0;l<szp;l++) for(m=0;m<szp;m++) if(pdf[k*szp*szp+l*szp+m]>0.0) cnt+=1.0;
  } else if(dim==4) {
    for(k=0;k<szp;k++) for(l=0;l<szp;l++) for(m=0;m<szp;m++) for(n=0;n<szp;n++) {
      if(pdf[k*szp*szp*szp+l*szp*szp+m*szp+n]>0.0) cnt+=1.0; }
  } else {
    printf("pdfnz WARNA: invalid dim=%d for pdf!\n",dim);
  }
  return cnt;
}








static double tetrospks3 (double* X1d,double* X2d,double* X3d,double* XO,int szX1,int szXO,int shuf) {
 
  double *pX3FX3P , *pX3P, *pX3FX1PX2PX3P, *pX1PX2PX3P, te, ds, tmp[5], mout[200], mean, norm, teout;
  double cnt1,cnt2,cnt3,jmax,N; 
  int i,j,k,l,m,n,sz,szp,*X1,*X2,*X3,minv1,maxv1,minv2,maxv2,minv3,maxv3;

  if (shuf>200) {printf("tetrospks3 INTERR nshuf (%d) >200\n",shuf); hxe();}
  if(useslice) { 
    if (end>0.0 && end<=szX1) szX1=(int)end; 
    printf("WARNING: using newly modified useslice capability\n");
  } else end=beg=0;
  sz=szX1-(int)beg; 
  X1=iscrset(sz*3); X2=X1+sz; X3=X1+2*sz;
  if(verbose>3) printf("X1:%p , X2:%p, X3:%p, scr:%p\n",X1,X2,X3,iscr);
  minv1=minv2=minv3=INT_MAX; maxv1=maxv2=maxv3=INT_MIN; cnt1=cnt2=cnt3=0;
  for (i=0;i<sz;i++) { 
    X1[i]=(int)X1d[i+(int)beg]; X2[i]=(int)X2d[i+(int)beg]; X3[i]=(int)X3d[i+(int)beg];
    if (X1[i]>0) cnt1++; if (X2[i]>0) cnt2++; if(X3[i]>0) cnt3++;
    if (X1[i]>maxv1) maxv1=X1[i];  if (X1[i]<minv1) minv1=X1[i]; 
    if (X2[i]>maxv2) maxv2=X2[i];  if (X2[i]<minv2) minv2=X2[i]; 
    if (X3[i]>maxv3) maxv3=X3[i];  if (X3[i]<minv3) minv3=X3[i]; 
  }
  if (minv1<0 || minv2<0 || minv3<0) {
    printf("tetrospks3 ERRB: minimum value must be >= 0:%d,%d,%d\n",minv1,minv2,minv3);hxe();}
  count+=1;
  if (minv1==maxv1 || minv2==maxv2 || minv3==maxv3) { 
    if(verbose>0) printf("tetrospks3 WARN0: #1:%d,%d,#2:%d,%d,#3:%d,%d)\n",minv1,maxv1,minv2,maxv2,minv3,maxv3);
    for (i=0;i<szXO;i++) XO[i]=0.0; if(szXO>=4+shuf)XO[shuf+3]=1.0; return 0.; }
  szp=(maxv1>maxv2)?(maxv1+1):(maxv2+1); if(maxv3+1>szp) szp=maxv3+1;
  if(verbose>1){printf("minv1:%d,maxv1:%d,cnt1:%g\n",minv1,maxv1,cnt1);
                printf("minv2:%d,maxv2:%d,cnt2:%g\n",minv2,maxv2,cnt2);
                printf("minv3:%d,maxv3:%d,cnt3:%g\n",minv3,maxv3,cnt3);}
  pX3P = (double*) calloc(szp,sizeof(double));
  pX3FX3P = (double*) calloc(szp*szp,sizeof(double));
  pX1PX2PX3P = (double*) calloc(szp*szp*szp,sizeof(double));
  pX3FX1PX2PX3P = (double*) calloc(szp*szp*szp*szp,sizeof(double));

  
  for(k=1;k<sz;k++) pX3P[ X3[k-1] ]++; 
  if(KTProb) { jmax=pdfnz(pX3P,szp,1); for(k=minv3;k<=maxv3;k++) if(pX3P[k]>0.0) {
      pX3P[k] = (0.5+pX3P[k]) / ( sz-1.0 + 0.5*jmax );}
  } else for(k=minv3;k<=maxv3;k++) pX3P[k] /= (sz-1);

  if (verbose>2) pdfpr(pX3P,szp,1,"pX3P");

  for(k=0;k<sz-1;k++) pX3FX3P[ X3[k+1]*szp + X3[k] ]++; 

  if(KTProb) { jmax=pdfnz(pX3FX3P,szp,2); 
    for(k=minv3;k<=maxv3;k++) for(l=minv3;l<=maxv3;l++) if(pX3FX3P[k*szp+l]>0.0){
      pX3FX3P[k*szp+l] = (pX3FX3P[k*szp+l]+0.5) / ( sz-1.0 + 0.5*jmax ); }
  } else for(k=minv3;k<=maxv3;k++) for(l=minv3;l<=maxv3;l++) pX3FX3P[k*szp+l]/=(sz-1);
  if (verbose>3) pdfpr(pX3FX3P,szp,2,"pX3FX3P");

  

  norm=entropxfgxpd(pX3P,pX3FX3P,minv3,maxv3,szp);               

  if (verbose>2) printf("H(X3F|X3P)=%g\n",norm);

  for (j=0,mean=0.;j<=shuf;j++) {

    
    if (j>0) { ishuffle(X1,sz); ishuffle(X2,sz); 
      memset(pX1PX2PX3P,0,sizeof(double)*szp*szp*szp); 
      memset(pX3FX1PX2PX3P,0,sizeof(double)*szp*szp*szp*szp); }
    for(l=1;l<sz;l++) pX1PX2PX3P[ X1[l-1]*szp*szp +X2[l-1]*szp + X3[l-1] ]++;  
    if(KTProb) {jmax=pdfnz(pX1PX2PX3P,szp,3); 
      for(l=minv1;l<=maxv1;l++) for(m=minv2;m<=maxv2;m++) for(n=minv3;n<=maxv3;n++) if(pX1PX2PX3P[l*szp*szp+m*szp+n]>0.0){
        pX1PX2PX3P[l*szp*szp+m*szp+n] = ( pX1PX2PX3P[l*szp*szp+m*szp+n] + 0.5 ) / ( sz-1.0 + 0.5*jmax ); }
    } else for(l=minv1;l<=maxv1;l++) for(m=minv2;m<=maxv2;m++) for(n=minv3;n<=maxv3;n++) pX1PX2PX3P[l*szp*szp+m*szp+n]/=(sz-1);
    if (verbose>3) pdfpr(pX1PX2PX3P,szp,3,"pX1PX2PX3P");

    
    for(k=0;k<sz-1;k++) pX3FX1PX2PX3P[ X3[k+1]*szp*szp*szp + X1[k]*szp*szp + X2[k]*szp + X3[k] ]++; 
    if(KTProb){jmax=pdfnz(pX3FX1PX2PX3P,szp,4);
      for(k=minv3;k<=maxv3;k++)for(l=minv1;l<=maxv1;l++)for(m=minv2;m<=maxv2;m++)for(n=minv3;n<=maxv3;n++){
        if(pX3FX1PX2PX3P[k*szp*szp*szp+l*szp*szp+m*szp+n]>0.0){
          pX3FX1PX2PX3P[k*szp*szp*szp+l*szp*szp+m*szp+n]=(pX3FX1PX2PX3P[k*szp*szp*szp+l*szp*szp+m*szp+n]+0.5)/(sz-1.0 + 0.5*jmax);}}
    } else for(k=minv3;k<=maxv3;k++) for(l=minv1;l<=maxv1;l++) for(m=minv2;m<=maxv2;m++) for(n=minv3;n<=maxv3;n++) {
      pX3FX1PX2PX3P[k*szp*szp*szp+l*szp*szp+m*szp+n]/=(sz-1);  }
    if (verbose>3) pdfpr(pX3FX1PX2PX3P,szp,4,"pX3FX1PX2PX3P");

    
    te = norm - entropx3fgx1px2px3pd(pX3FX1PX2PX3P,pX1PX2PX3P,minv1,maxv1,minv2,maxv2,minv3,maxv3,szp);
    if (j>0) {mean+=te; mout[j-1]=te;} else teout=te; 
  }
  if (shuf>0) mean/=shuf;
  if (szXO>0) XO[0]=teout; if (szXO>1) XO[1]=norm; if (szXO>2) XO[2]=mean;
  if (szXO>=3+shuf) for (i=0;i<shuf;i++) XO[i+3]=mout[i];
  if (szXO>=4+shuf && shuf>0) { 
    cnt1 = 0.0; 
    for(i=0;i<shuf;i++) if(teout < XO[i+3]) {cnt1+=1.0; if(verbose>2)printf("teout:%g, XO[%d]=%g\n",teout,i+3,XO[i+3]);}
    XO[shuf+3] = cnt1 / (double)shuf;
    if(verbose>2) printf("cnt1=%g, shuf=%d, XO[%d]=%g\n",cnt1,shuf,shuf+3,XO[shuf+3]);
  }
  if(verbose>2) printf("te=%g\n",te);
  free(pX3P); free(pX3FX3P); free(pX1PX2PX3P); free(pX3FX1PX2PX3P);
  return (teout-mean)/norm;
}



    
    
    
static double tetrospks2 (double* X1d,double* X2d,double* XO,int szX1,int szXO,int shuf) {
  double *pX2P,*pX2FX2PX1P,*pX2PX1P,*pX2FX2P,te,ds,tmp[5],mout[200],mean,norm,teout;
  double cnt1,cnt2,jmax,N; int i,j,k,l,m,sz,szp,*X1,*X2,minv1,maxv1,minv2,maxv2;
  if (shuf>200) {printf("tetrospks2 INTERR nshuf (%d) >200\n",shuf); hxe();}
  if(useslice) { 
    if (end>0.0 && end<=szX1) szX1=(int)end; 
    printf("WARNING: using newly modified useslice capability\n");
  } else end=beg=0;
  sz=szX1-(int)beg; 
  X1=iscrset(sz*2); X2=X1+sz; 
  if(verbose>3) printf("X1:%p , X2:%p, scr:%p\n",X1,X2,iscr);
  minv1=minv2=INT_MAX; maxv1=maxv2=INT_MIN; cnt1=cnt2=0;
  if(verbose>2) printf("before: minv1=%d ,maxv1=%d, minv2=%d, maxv2=%d\n",minv1,maxv1,minv2,maxv2);
  for (i=0;i<sz;i++) { 
    X1[i]=(int)X1d[i+(int)beg]; X2[i]=(int)X2d[i+(int)beg];
    if (X1[i]>0) cnt1++;           if (X2[i]>0) cnt2++; 
    if (X1[i]>maxv1) maxv1=X1[i];  if (X1[i]<minv1) minv1=X1[i]; 
    if (X2[i]>maxv2) maxv2=X2[i];  if (X2[i]<minv2) minv2=X2[i]; 
  }
  if(verbose>2) printf("after: minv1=%d ,maxv1=%d, minv2=%d, maxv2=%d\n",minv1,maxv1,minv2,maxv2);
  if (minv1<0 || minv2<0) {
    printf("tentropspks ERRB: minimum value must be >= 0:%d,%d\n",minv1,minv2);hxe();}
  if (binmin) { 
    cnt1/=sz; cnt2/=sz; 
    if (cnt1<binmin) return -11.; else if (cnt2<binmin) return -12.;
    if (binmax) if (cnt1>binmax) return -11.; else if (cnt2>binmax) return -12.;
    if (abs(cnt1-cnt2)>cutoff) return -13.;
  }
  if(verbose>2)printf("tentropspks:minv1=%d,maxv1=%d,minv2=%d,maxv2=%d\n",minv1,maxv1,minv2,maxv2);
  count+=1;
  if (minv1==maxv1 || minv2==maxv2) { 
    if(verbose>0) printf("tentropspk WARN0: #1:%d,%d,#2:%d,%d)\n",minv1,maxv1,minv2,maxv2);
    for (i=0;i<szXO;i++) XO[i]=0.0; if(szXO>=4+shuf)XO[shuf+3]=1.0; return 0.; }
  szp=(maxv1>maxv2)?(maxv1+1):(maxv2+1);
  pX2P = (double*) calloc(szp,sizeof(double));
  pX2PX1P = (double*) calloc(szp*szp,sizeof(double));
  pX2FX2P = (double*) calloc(szp*szp,sizeof(double));
  pX2FX2PX1P = (double*) calloc(szp*szp*szp,sizeof(double));
  
  for(k=1;k<sz;k++) pX2P[ X2[k-1] ]++; 
  if(KTProb) { jmax=pdfnz(pX2P,szp,1); for(k=minv2;k<=maxv2;k++) if(pX2P[k]>0.0) {
      pX2P[k] = (0.5+pX2P[k]) / ( sz-1.0 + 0.5*jmax );}
  } else for(k=minv2;k<=maxv2;k++) pX2P[k] /= (sz-1);
  if (verbose>2) pdfpr(pX2P,szp,1,"pX2P");
  for(k=0;k<sz-1;k++) pX2FX2P[ X2[k+1]*szp + X2[k] ]++; 
  if(KTProb) { jmax=pdfnz(pX2FX2P,szp,2); 
    for(k=minv2;k<=maxv2;k++) for(l=minv2;l<=maxv2;l++) if(pX2FX2P[k*szp+l]>0.0){
      pX2FX2P[k*szp+l] = (pX2FX2P[k*szp+l]+0.5) / ( sz-1.0 + 0.5*jmax ); }
  } else for(k=minv2;k<=maxv2;k++) for(l=minv2;l<=maxv2;l++) pX2FX2P[k*szp+l]/=(sz-1);
  if (verbose>3) pdfpr(pX2FX2P,szp,2,"pX2FX2P");
  
  norm=entropxfgxpd(pX2P,pX2FX2P,minv2,maxv2,szp);               
  if (verbose>2) printf("H(X2F|X2P)=%g\n",tmp[0]);
  for (j=0,mean=0.;j<=shuf;j++) {
    
    if (j>0) { ishuffle(X1,sz); 
      memset(pX2PX1P,0,sizeof(double)*szp*szp); 
      memset(pX2FX2PX1P,0,sizeof(double)*szp*szp*szp);    }
    for(l=1;l<sz;l++) pX2PX1P[ X2[l-1]*szp + X1[l-1] ]++;  
    if(KTProb) {jmax=pdfnz(pX2PX1P,szp,2); 
      for(l=minv2;l<=maxv2;l++) for(m=minv1;m<=maxv1;m++) if(pX2PX1P[l*szp+m]>0.0){
        pX2PX1P[l*szp+m] = ( pX2PX1P[l*szp+m] + 0.5 ) / ( sz-1.0 + 0.5*jmax ); }
    } else for(l=minv2;l<=maxv2;l++) for(m=minv1;m<=maxv1;m++) pX2PX1P[l*szp+m]/=(sz-1);
    if (verbose>3) pdfpr(pX2PX1P,szp,2,"pX2PX1P");
    
    for(k=0;k<sz-1;k++) pX2FX2PX1P[ X2[k+1]*szp*szp + X2[k]*szp + X1[k] ]++; 
    if(KTProb){jmax=pdfnz(pX2FX2PX1P,szp,3);
      for(k=minv2;k<=maxv2;k++)for(l=minv2;l<=maxv2;l++)for(m=minv1;m<=maxv1;m++){
        if(pX2FX2PX1P[k*szp*szp+l*szp+m]>0.0){
          pX2FX2PX1P[k*szp*szp+l*szp+m]=(pX2FX2PX1P[k*szp*szp+l*szp+m]+0.5)/(sz-1.0 + 0.5*jmax);}}
    } else for(k=minv2;k<=maxv2;k++) for(l=minv2;l<=maxv2;l++) for(m=minv1;m<=maxv1;m++) {
      pX2FX2PX1P[k*szp*szp+l*szp+m]/=(sz-1);  }
    if (verbose>3) pdfpr(pX2FX2PX1P,szp,3,"pX2FX2PX1P");
    
    te = norm - entropx2fgx2px1pd(pX2FX2PX1P,pX2PX1P,minv1,maxv1,minv2,maxv2,szp);
    if (j>0) {mean+=te; mout[j-1]=te;} else teout=te; 
  }
  if (shuf>0) mean/=shuf;
  if (szXO>0) XO[0]=teout; if (szXO>1) XO[1]=norm; if (szXO>2) XO[2]=mean;
  if (szXO>=3+shuf) for (i=0;i<shuf;i++) XO[i+3]=mout[i];
  if (szXO>=4+shuf && shuf>0) { 
    cnt1 = 0.0; 
    for(i=0;i<shuf;i++) if(teout < XO[i+3]) {cnt1+=1.0; if(verbose>2)printf("teout:%g, XO[%d]=%g\n",teout,i+3,XO[i+3]);}
    XO[shuf+3] = cnt1 / (double)shuf;
    if(verbose>2) printf("cnt1=%g, shuf=%d, XO[%d]=%g\n",cnt1,shuf,shuf+3,XO[shuf+3]);
  }
  if(verbose>2) printf("te=%g\n",te);
  free(pX2FX2P); free(pX2PX1P); free(pX2P); free(pX2FX2PX1P);
  return (teout-mean)/norm;
}


static double pdfpr (double* pdf,int szp,int dim, char* name) {
  double x,ds; int i,j,k,l,m,cnt,*nonzero;
  ds=0.; 
  printf("Contents of PDF %s\n",name);
  if (dim>2) { 
    nonzero=(int *)calloc(szp,sizeof(int));
    for(k=0,cnt=0;k<szp;k++) { 
      for(l=0;l<szp;l++) for (m=0;m<szp;m++) if (pdf[k*szp*szp+l*szp+m]>0.) cnt++;
      if (cnt>0) nonzero[k]=1; 
    }
  }
  if (dim==1) {
    for(m=0;m<szp;m++) { printf("%g ",pdf[m]); ds+=pdf[m]; }
  } else if (dim==2) {
    for(l=0;l<szp;l++) { 
      for(m=0;m<szp;m++) { printf("%g ",pdf[l*szp+m]); ds+=pdf[l*szp+m]; }
      printf("\n");
    }
  } else { 
    for (k=0;k<szp;k++) if (nonzero[k]==1) { 
      printf("\tSlice #%d:\n",k); 
      for(l=0;l<szp;l++) { 
        for(m=0;m<szp;m++) { printf("%g ",(x=pdf[k*szp*szp+l*szp+m])); ds+=x; }
        printf("\n");
      }
    }
  }
  printf("ds (%s) = %g\n",name,ds);
}



static double getbandwidth (void* vv) {
  double *x, sd;
  int n;
  n = vector_instance_px(vv,&x);
  sd = getstdevd(x,n);
  return getbandwidthd(1,n,sd);
}



static double getdisc (void* vv) {
  double *pin,*pout;
  int szin,szout,nbins,*ptmp,i;
  szin = vector_instance_px(vv,&pin);
  if((szout = vector_arg_px(1,&pout))<szin) {
    printf("getnorm ERRA: output size %d < input size %d\n",szout,szin);
    return 0.0;
  }
  nbins = (int)*getarg(2);
  if(nbins < 1) {
    printf("getnorm ERRB: nbins must be positive\n");
    return 0.0;
  }
  ptmp = getnormd(pin,szin,nbins);
  for(i=0;i<szin;i++) pout[i] = ptmp[i];
  free(ptmp);
  return 1.0;
}



static double downsampavg (void* vv) {
  double *pin,*pout,*ptmp,dsz;
  int szin,szout,winsz;
  szin=vector_instance_px(vv,&pin);
  winsz = (int)*getarg(2);
  if(winsz < 1) {
    printf("downsampavg ERRA: winsz must be positive\n");
    return 0.0;
  }
  dsz= (double)szin/(double)winsz;  dsz = ceil(dsz);
  if((szout=vector_arg_px(1,&pout))<dsz) {
    printf("downsampavg ERRB: output size %d < required size %d\n",szout,(int)dsz);
    return 0.0;
  }
  if(!downsampavgd(pin,pout,szin,szout,winsz)) {
    printf("downsampavg ERRC: output size too small\n");
    return 0.0;
  }
  vector_resize(vector_arg(1),(int)dsz);
  return 1.0;
}


static double chgenc (void* vv) {
  double *pin,*pout;
  int szin,szout;
  szin=vector_instance_px(vv,&pin);
  if((szout=vector_arg_px(1,&pout))!=szin) {
    printf("chgenc ERRA: output size %d < required size %d\n",szout,szin);
    return 0.0;
  }
  chgencd(pin,pout,szin);
  return 1.0;
}



static double kprob1D (void* vv) {
  double *x,val,h;
  int n;
  n = vector_instance_px(vv,&x);  h = *getarg(1);  val = *getarg(2);
  return kprob1Dd(x,n,h,val);
}



static int ITcompare(const void* a, const void* b)
{ const int i1 = *(const int*)a;
  const int i2 = *(const int*)b;
  const double term1 = ITsortdata[i1];
  const double term2 = ITsortdata[i2];
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}


void ITcsort (int n, const double mdata[], int index[])

{ int i;
  ITsortdata = mdata;
  for (i = 0; i < n; i++) index[i] = i;
  qsort(index, n, sizeof(int), ITcompare);
}


static double* ITgetrank (int n, double mdata[])





{ int i;
  double* rank;
  int* index;
  rank = calloc(n,sizeof(double));
  if (!rank) return NULL;
  index = calloc(n,sizeof(int));
  if (!index)
  { free(rank);
    return NULL;
  }
  
  ITcsort (n, mdata, index);
  
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



static double mutinfbd(double* x,double* y,int sz,int nbins) {
  double *xr,*yr,*px,*py,**pxy,ret,dinf;
  int i,j,idx,idy;
  ret=-1.0;
  xr=yr=px=py=NULL; pxy=NULL;
  if(verbose>0) printf("sz=%d, nbins=%d\n",sz,nbins);
  if(!(xr=ITgetrank(sz,x))) { printf("mutinfb ERRB: out of memory\n"); goto MUTEND; }
  if(!(yr=ITgetrank(sz,y))) { printf("mutinfb ERRB: out of memory\n"); goto MUTEND; }
  if(verbose>0) { printf("xr: "); for(i=0;i<sz;i++) printf("%g ",xr[i]); printf("\n");
                  printf("yr: "); for(i=0;i<sz;i++) printf("%g ",yr[i]); printf("\n");}
  if(!(px=(double*)calloc(nbins,sizeof(double)))){printf("mutinfb ERRB: out of memory\n"); goto MUTEND;}
  if(!(py=(double*)calloc(nbins,sizeof(double)))){printf("mutinfb ERRB: out of memory\n"); goto MUTEND;}
  if(!(pxy=getdouble2D(nbins,nbins))){printf("mutinfb ERRB: out of memory\n"); goto MUTEND;}
  for(i=0;i<sz;i++) { 
    idy = (int) nbins * ( yr[i] / (sz-1) ); if(idy>=nbins) idy = nbins-1;
    idx = (int) nbins * ( xr[i] / (sz-1) ); if(idx>=nbins) idx = nbins-1;
    px[idx]++; py[idy]++; pxy[idy][idx]++;
  }
  if(verbose>0) { 
    printf("px: "); for(i=0;i<nbins;i++) printf("%.2f ",px[i]); printf("\n");
    printf("py: "); for(i=0;i<nbins;i++) printf("%.2f ",py[i]); printf("\n");
    printf("pxy: \n");
    for(i=0;i<nbins;i++) { for(j=0;j<nbins;j++) printf("%.2f ",pxy[i][j]); printf("\n"); }
  }
  dinf=0.0;
  for(idy=0;idy<nbins;idy++) { 
    if(py[idy]==0.0) continue;
    for(idx=0;idx<nbins;idx++) {
      if(px[idx]==0.0 || pxy[idy][idx]==0.0) continue;
      dinf += pxy[idy][idx] * log2d( pxy[idy][idx] / ( px[idx] * py[idy] ) );
    }
  }
  if(verbose>0) printf("A: dinf = %g\n",dinf);
  dinf = dinf / (double) sz;
  if(verbose>0) printf("B: dinf = %g\n",dinf);
  dinf += log2d(sz);
  if(verbose>0) printf("C: dinf = %g\n",dinf);
  ret = dinf;
MUTEND:
  if(xr) free(xr); if(yr) free(yr);  if(px) free(px); if(py) free(py);
  if(pxy) freedouble2D(&pxy,nbins);
  return ret;
}


static double mutinfb (void* v) {
  double *x,*y;
  int sz,nbins,tmp;
  sz = vector_instance_px(v,&x);
  if((tmp=vector_arg_px(1,&y))!=sz) {
    printf("mutinfb ERRA: x,y must have same sizes : %d %d\n",sz,tmp);
    return -1.0;
  }
  nbins = ifarg(2) ? (int) *getarg(2) : 10;
  return mutinfbd(x,y,sz,nbins);
}

double entropspksd (double* x,int sz) {
  double *px,ret,dinf,size;
  int i,*X,maxv,minv,cnt,err;
  if(sz<1) {printf("entropspks ERR0: min size must be > 0!\n"); return -1.0;}
  size=(double)sz;
  ret=-1.0; err=0;
  px=NULL; 
  minv=1000000000; maxv=-1000000000;
  X=scrset(sz*2); 
  cnt=0;
  for (i=0;i<sz;i++) { 
    X[i]=(int)x[i]; 
    if (X[i]>0) cnt++;  
    if (X[i]>maxv) maxv=X[i];  if (X[i]<minv) minv=X[i]; 
  }
  if (minv<0){printf("entropb ERRB: minimum value must be >= 0:%d\n",minv);hxe();}
  if(!(px=(double*)calloc(maxv+1,sizeof(double)))){
    printf("entropb ERRB: out of memory\n"); err=1; goto ENTEND;}
  for(i=0;i<sz;i++) px[X[i]] += 1.0; 
  for(i=minv;i<=maxv;i++) px[i]/=size;
  if(verbose>1) {printf("px: "); for(i=minv;i<=maxv;i++) printf("%.2f ",px[i]); printf("\n");}
  dinf=0.0;
  if(usetable) {
    for(i=minv;i<=maxv;i++) if(px[i]>0.) {
      if(px[i]>=MINLOG2&&px[i]<=MAXLOG2) dinf -= px[i]*_n_LOG2(px[i]); else {
                                         dinf -= px[i]*  log2d(px[i]); }
    }
  } else for(i=minv;i<=maxv;i++) if(px[i]>0.) dinf -= px[i] * log2d(px[i]);
  ret = dinf;
ENTEND:
  if(px) free(px);
  if (err) hxe();
  return ret;
}

static double entropspks (void* v) {
  double *x;
  int sz;
  sz = vector_instance_px(v,&x);
  return entropspksd(x,sz);
}

static double mxentropd (double* m, int rows, int cols, int tsz) {
  int lsz,i,j,*pcnts,x,y,tdx,yy,xx,ntiles,*ptmp;   double dent,p,d;
  lsz=1; ntiles=0; dent=0.0;
  for(i=0;i<tsz*tsz;i++) lsz*=2;
  if(verbose>=1) printf("# of tile patterns: %d\n",lsz);
  ptmp=verbose>=1?(int*)calloc(sizeof(int),tsz*tsz):0x0;
  pcnts = (int*) calloc(sizeof(int),lsz);
  for(y=0;y<rows-tsz;y++) {
    for(x=0;x<cols-tsz;x++) {
      tdx=0; i=1; j=0;
      for(yy=y;yy<y+tsz;yy++) {
        for(xx=x;xx<x+tsz;xx++) {
          d = m[yy*cols+xx]; 
          if(ptmp) ptmp[j]=(int)d; 
          if(d) tdx += i;
          i *= 2; j+=1;
        }
      }
      pcnts[tdx]++; ntiles+=1;
      if(ptmp) {
        for(j=0;j<tsz*tsz;j++) printf("%d",ptmp[j]);
        printf("\tidx:%d\n",tdx); 
      }
    }
  }
  for(i=0;i<lsz;i++) {
    if(verbose>=1) printf("pcnts[%d]=%d\n",i,pcnts[i]);
    if(pcnts[i]) {
      p = (double) pcnts[i] / (double) ntiles;
      dent -= p*log2d(p);
    }
  }
  free(pcnts); if(ptmp) free(ptmp);
  return dent;
}




static double mxentrop (void* v) {
  double* m;
  int sz,rows,cols,tsz;
  sz = vector_instance_px(v,&m);
  rows = (int)*getarg(1);
  cols = sz/rows;
  tsz = (int)*getarg(2);
  return mxentropd(m,rows,cols,tsz);
}

ENDVERBATIM


FUNCTION testoddometer () {
  VERBATIM
  int *po,i;
  po = (int*) calloc(4,sizeof(int));
  do {
    for(i=0;i<4;i++) printf("%d",po[i]);
    printf("\n");
  } while(oddometerinc(po,4,9));
  free(po);
  return 1.0;
  ENDVERBATIM
}

FUNCTION EXP (x) {
 TABLE DEPEND MINEXP,MAXEXP FROM MINEXP TO MAXEXP WITH 50000
 EXP = exp(x)
}

FUNCTION LOG2 (x) {
 TABLE DEPEND MINLOG2,MAXLOG2 FROM MINLOG2 TO MAXLOG2 WITH 50000
 LOG2 = log(x)/LG2
}

PROCEDURE install () {
  if (installed==1) {
    printf("infot.mod version %s\n","$Id
  } else {
    installed = 1
    VERBATIM
    _check_LOG2(); _check_EXP();
    install_vector_method("tentrop",tentrop);
    install_vector_method("ntedir",ntedir);
    install_vector_method("getbandwidth",getbandwidth);
    install_vector_method("getdisc",getdisc);
    install_vector_method("downsampavg",downsampavg);
    install_vector_method("chgenc",chgenc);
    install_vector_method("kprob1D",kprob1D);
    install_vector_method("mutinfb",mutinfb);
    install_vector_method("tentropspks",tentropspks);
    install_vector_method("nte",nte);
    install_vector_method("ntel2",ntel2);
    install_vector_method("entropspks",entropspks);
    install_vector_method("entropxfgxp",entropxfgxp);
    install_vector_method("lz76c",lz76c);
    install_vector_method("mxentrop",mxentrop);
    ENDVERBATIM
  }
}