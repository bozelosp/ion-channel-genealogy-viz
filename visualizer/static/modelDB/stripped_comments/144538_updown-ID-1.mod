NEURON {
  SUFFIX nothing
  
  GLOBAL UPDOWN_INSTALLED, SHM_UPDOWN, NOV_UPDOWN, DEBUG_UPDOWN
}

PARAMETER {
  UPDOWN_INSTALLED=0
  DEBUG_UPDOWN=0
  SHM_UPDOWN=4   
  NOV_UPDOWN=1   
  CREEP_UPDOWN=0 
}

VERBATIM
#include "misc.h"
static void hxf(void *ptr) { free(ptr); hoc_execerror("",0); }
ENDVERBATIM







VERBATIM

  
#define UDSL 500
#define UDNQ 11

#define LOC     nq[0] 
#define PEAK  	nq[1] 
#define WIDTH  	nq[2] 
#define BASE  	nq[3] 
#define HEIGHT  nq[4] 
#define START  	nq[5] 
#define SLICES  nq[6] 
#define SHARP  	nq[7] 
#define INDEX  	nq[8] 

#define NESTED  nq[10] 
  
static double updown (void* vv) {
  int i, k, m, n, nqsz, nsrc, jj[UDSL], f[UDSL], lc, dsz[UDSL], nqmax, thsz, lc2, done, dbn;
  double *src, *tvec, *th, *dest[UDSL], *nq[UDNQ], *tmp, *dbx, lt, thdist;
  Object *ob, *ob2;
  void *vvd[UDSL], *vvth, *vnq[UDNQ];
  
  nsrc = vector_instance_px(vv, &src); 
  thsz = vector_arg_px(1, &th);        
  ob =  *hoc_objgetarg(2);             
  ob2 = *hoc_objgetarg(3);             
  tmp = (double *)ecalloc(nsrc, sizeof(double));  
  lc =  ivoc_list_count(ob);
  lc2 = ivoc_list_count(ob2);
  if (lc>UDSL) {printf("updown ERRF mismatch: max slice list:%d %d\n",UDSL,lc); hxf(tmp);}
  if (lc2!=UDNQ){printf("updown ERRB mismatch: NQS sz is %d (%d in list)\n",UDNQ,lc2);hxf(tmp);}
  if (nsrc<lc) {printf("updown ERRC mismatch: %d %d\n",lc,nsrc); hxf(tmp);} 
  if (lc!=thsz) {printf("updown ERRA mismatch: %d %d\n",lc,thsz); hxf(tmp);}
  if (!ismono1(th,thsz,-1)) {printf("updown ERRD: not mono dec %g %d\n",th[0],thsz); hxf(tmp);}
  
  for (k=0;k <lc;k++)  dsz[k] =list_vector_px3(ob , k, &dest[k], &vvd[k]);
  for (k=0;k<lc2;k++) {
    i=list_vector_px3(ob2, k, &nq[k],   &vnq[k]);
    if (k==0) nqmax=i; else if (i!=nqmax) { 
      printf("updown ERRE mismatch: %d %d %d\n",k,i,nqmax); hxf(tmp); }
  }
  
  
  
  for (k=0; k<lc; k++) {   
    jj[k]=f[k]=0; 
    for (i=0;i<nsrc && src[i]>th[k];i++) {} 
    for (; i<nsrc; i++) { 
      if (src[i]>th[k]) { 
        if (f[k]==0) { 
          if (jj[k]>=dsz[k]){printf("(%d,%d,%d) :: ",k,jj[k],dsz[k]);
            hoc_execerror("Dest vec too small in updown ", 0); }
          dest[k][jj[k]++] = (i-1) + (th[k]-src[i-1])/(src[i]-src[i-1]); 
          f[k]=1; 
          tmp[k]=-1e9; dest[k][jj[k]]=-1.; 
        }
        if (f[k]==1 && src[i]>tmp[k]) { 
          tmp[k]=src[i]; 
          dest[k][jj[k]] = (double)i; 
        }
      } else {          
        if (f[k]==1) {  
          jj[k]++;      
          dest[k][jj[k]++] = (i-1) + (src[i-1]-th[k])/(src[i-1]-src[i]);
          f[k]=0; 
        }
      }
    }
  }
  
  for (k=0;k<lc;k++) vector_resize(vvd[k],(int)(floor((double)jj[k]/3.)*3.));
  for (i=0; i<nsrc; i++) tmp[i]=0.; 
  
  
  
  for (k=0;k<lc;k++) { 
    for (i=1;i<jj[k];i+=3) { 
      m=(int)dest[k][i]; 
      if (tmp[m-2]<0 || tmp[m-1]<0 || tmp[m+1]<0 || tmp[m+2]<0) continue; 
      tmp[m]--;  
      tmp[m-1]=dest[k][i-1]; tmp[m+1]=dest[k][i+1]; 
    }
  }
  
  
  
  
  for (i=0,k=0; i<nsrc; i++) if (tmp[i]<0.) { 
    if (k>=nqmax) { printf("updown ERRG OOR in NQ db: %d %d\n",k,nqmax); hxf(tmp); }
    LOC[k]=(double)i;  
    WIDTH[k]=tmp[i+1]; 
    START[k]=tmp[i-1]; 
    SLICES[k]=-tmp[i];  
    k++;
  }
  nqsz=k;   
  if (DEBUG_UPDOWN && ifarg(4)) { dbn=vector_arg_px(4, &dbx); 
    if (dbn<nsrc) printf("updown ERRH: Insufficient room in debug vec (%d<%d)\n",dbn,nsrc); 
    else for (i=0;i<nsrc;i++) dbx[i]=tmp[i]; 
  }
  
  
  
  
  
  
  
  
  
  
  
  if (NOV_UPDOWN) for (i=0;i<nqsz;i++) { 
    if ((i-1)>0 && START[i] < LOC[i-1]) { 
      if (DEBUG_UPDOWN) printf("LT problem %d %g %g<%g\n",i,LOC[i],START[i],LOC[i-1]);
      for (m=lc-1,done=0;m>=0 && !done;m--) { 
        for (n=1;n<jj[m] && !done;n+=3) {     
          
          if (floor(dest[m][n])==LOC[i] && dest[m][n-1]>LOC[i-1]) {
            
            
            START[i]=dest[m][n-1]; WIDTH[i]=dest[m][n+1]; done=1; 
          }
        }
      }
    }
    
    if ((i+1)<nqsz && WIDTH[i]>LOC[i+1]) {
      if (DEBUG_UPDOWN) printf("RT problem %d %g %g>%g\n",i,LOC[i],WIDTH[i],LOC[i+1]);
      for (m=lc-1,done=0;m>=0 && !done;m--) { 
        for (n=1;n<jj[m] && !done;n+=3) {     
          
          if (floor(dest[m][n])==LOC[i] && dest[m][n+1]<LOC[i+1]) {
            
            START[i]=dest[m][n-1]; WIDTH[i]=dest[m][n+1]; done=1;
          }
        }
      }        
    }
  }

  
  
  
  
  if(CREEP_UPDOWN) for(i=0,k=0;i<nsrc;i++) if(tmp[i]<0.){

    
    int idx = (int)START[k];
    while(idx >= 1 && src[idx] >= src[idx-1]) idx--;
    START[k] = idx;

    
    idx = (int)WIDTH[k];
    while(idx < nsrc-1 && src[idx] >= src[idx+1]) idx++;
    WIDTH[k] = idx;

    k++;
  }

  
  
  
  for (i=0,k=0; i<nsrc; i++) if (tmp[i]<0.) { 
    
    lt=src[(int)floor(START[k])]+(START[k]-floor(START[k]))*\
      (src[(int)floor(START[k]+1.)]-src[(int)floor(START[k])]);
    BASE[k]=lt;         
    PEAK[k]=src[i];     
    WIDTH[k] = WIDTH[k] - START[k]; 
    HEIGHT[k]=PEAK[k]-BASE[k]; 
    
    
    
    SHARP[k]=(src[i]-src[i-(int)SHM_UPDOWN])-(src[i+(int)SHM_UPDOWN]-src[i]);
    INDEX[k]=(double)k;
    k++;
  }
  
  int iNumBumps = k;

  
  if(!NOV_UPDOWN){
    for(i=0; i<iNumBumps; i++){
      NESTED[i] = 0;
      int j = 0;
      for(;j<iNumBumps;j++){
        if(i!=j && LOC[j] >= START[i] && LOC[j] <= START[i]+WIDTH[i]){
          NESTED[i]+=1.0;
        }
      }
    }
  } else for(i=0;i<iNumBumps;i++) NESTED[i]=0.0;

  
  for (i=0;i<lc2;i++) vector_resize(vnq[i], nqsz);
  if (k!=nqsz) { printf("updown ERRI INT ERR: %d %d\n",k,nqsz); hxf(tmp); }
  free(tmp);
  return jj[0];
}


ENDVERBATIM


PROCEDURE install_updown () {
  if (UPDOWN_INSTALLED==1) {
    printf("$Id
  } else {
  UPDOWN_INSTALLED=1
  VERBATIM {
  install_vector_method("updown", updown);
  }
  ENDVERBATIM
  }
}