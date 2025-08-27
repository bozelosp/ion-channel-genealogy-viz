NEURON {
  SUFFIX staley
  GLOBAL installed,verbose,samprate,mindist,minspikes,abovebth
  GLOBAL Lintdur,Bintdur,spkup
}


PARAMETER {
  installed=0
  verbose=0
  samprate=2000    
  Lintdur=120      
  Bintdur=12       
  Sintdur=2        
  mindist=0.1      
  minspikes=40     
  abovebth=0       
  endwting=1       
  flag=0           
  spkup=12         
  spklim=0.25      
  useavgtots=0     

  usesharp=0       
  sharpoff=4       
  sharpth=-500     

  incby1=0         
}

VERBATIM

static double gszspk (double* channelData, int* spks, double* tots, int channelsLength,
                      int intervalLength);
static double gsz (double* channelData,double* diffcor,double* percentc,double* tots,
                   int channelsLength,int intervalLength, int SgroupCount,int LintCnt);
#include "misc.h"
static  ListVec* pL;

typedef struct SEIZURE {
  int startIntervalIndex;
  int endIntervalIndex;
  int totalSpikesCount;
  int endSpikesCount;
  int ID;
  int startIndex;
  int endIndex;
} sSeizure;





sSeizure* AllocSeizure (int intervalIndex,int spikesCount)
{
  sSeizure* p;
  if(!(p = (sSeizure*)calloc(1,sizeof(sSeizure)))){
    printf("AllocSeizure ERR: out of memory!\n"); hxe();
  }
  p->startIntervalIndex = p->endIntervalIndex = intervalIndex;
  p->endIntervalIndex+=1; 
  p->totalSpikesCount = spikesCount;
  p->ID = -1; 
  return p;
}

typedef struct LSEIZURE {
  int bufsz;
  int count;
  sSeizure** pp;
} LSeizure;

int InitLSeizure (LSeizure* pl,int sz) {
  pl->bufsz = sz;
  pl->count = 0;
  if(sz==0) pl->pp=NULL; else pl->pp=(sSeizure**)calloc(sz,sizeof(sSeizure*));
  return pl->bufsz;
}

void FreeLSeizure(LSeizure* pl) {
  int i;
  
  free(pl->pp);
}

int AddSeizure(LSeizure* pl, sSeizure* ps) {
  if(0) printf("pl=%p , pl->count=%d, pl->bufsz=%d\n",pl,pl->count,pl->bufsz);
  if( pl->count + 1 >= pl->bufsz ) {
    pl->bufsz *= 4;
    if(0)printf("pl=%p, realloc pl->count=%d, pl->bufsz=%d\n",pl,pl->count,pl->bufsz);
    if(! (pl->pp=(sSeizure**)realloc(pl->pp,sizeof(sSeizure*)*pl->bufsz)) ) {
      printf("AddSeizure: out of memory!\n"); hxe();
    }
  }
  pl->pp[pl->count++] = ps;
  return 1;
}

int RemoveSeizureAt(LSeizure* pl, int idx) {
  int i,j;
  if( idx < 0 || idx >= pl->count) {
    printf("RemoveSeizureAt: invalid index=%d, count=%d!\n",idx,pl->count); hxe();
  }
  for(i=idx+1;i<pl->count;i++) pl->pp[i-1]=pl->pp[i];
  pl->count--;
  return 1;
}

void printSeizure (sSeizure* p) {
  printf("ID:%d, startIntervalIndex:%d, endIntervalIndex:%d, totalSpikes:%d, endSpikes:%d, startIndex:%d, endIndex:%d\n",
	 p->ID,p->startIntervalIndex,p->endIntervalIndex,p->totalSpikesCount,p->endSpikesCount,p->startIndex,p->endIndex);
}

void printSeizures (LSeizure* lp) {
  int i;
  for(i=0;i<lp->count;i++) printSeizure(lp->pp[i]);
}


static LSeizure* dgetseizures (double* channelData,int channelsLength) {
  
  
  
  
  
  int intervalIndex, dataIndex, seizureIndex, iS, stopIndex, uplim;
  int intervalLength, channelIndex, i,j, intervalsCount, SgroupCount, true_var;
  int seizuresCount, spikesCount, totsLen, bufszStart, dbxi, ii, jj, LintCnt;
  double value, min_minS, max_maxS, HV, LV; 
  double diffc, totsum, diffthresh;
  int upflag, upcount, avgnumspikes, *spks;
  double upvalue, limit, *diffcor, *percentc, *tots, sharp;
  LSeizure tmpSeizures, *pSeizuresOut;  
  sSeizure *tmpSeizure, *currentSeizure, *nextSeizure; 
  tmpSeizure=currentSeizure=nextSeizure=0x0;
  diffcor=percentc=tots=0x0; spks=0x0;
  if(!(pSeizuresOut=(LSeizure*)calloc(1,sizeof(LSeizure)))){printf("getseizures ERR: out of memory!\n");hxe();}
  bufszStart=400;
  InitLSeizure(pSeizuresOut,bufszStart); InitLSeizure(&tmpSeizures,bufszStart);
  intervalLength = (int)Bintdur*samprate;
  intervalsCount = channelsLength / intervalLength; 
  uplim=(int)(0.5+spkup*samprate/1e3); 
  if(intervalsCount<1) {
    printf("getseizures ERR:invalid intervalsCount:%d %d\n",channelsLength,intervalLength);hxe();}
  SgroupCount = (int)Bintdur/Lintdur*1e3; 
  LintCnt=(int)Lintdur*samprate/1000;
  if(intervalsCount<1 || SgroupCount<1) { 
    printf("Error: Data length too short, cannot process."); hxe(); }
  if((diffcor = (double*) calloc(intervalsCount,sizeof(double)))==0x0) {
    printf("getseizures ERRdfc: out of memory!\n"); hxe();  }
  if(!(percentc = (double*) calloc(intervalsCount,sizeof(double)))) { 
    printf("getseizures ERRpcc: out of memory!\n"); hxe();  }            
  if(!(tots = (double*) calloc(intervalsCount,sizeof(double)))) {
    printf("getseizures ERRtts: out of memory!\n"); hxe();  }
  if(!(spks=(int*) calloc((size_t)(jj=intervalsCount*Bintdur/Sintdur),sizeof(int)))) {
    printf("getseizures ERRspk: out of memory!\n"); hxe();  }
  for (ii=0;ii<jj;ii++) spks[ii]=0; 
  totsLen = intervalsCount;
  seizuresCount=0;
  
  diffthresh=\
    gsz(channelData,diffcor,percentc,tots, channelsLength, intervalLength,SgroupCount,LintCnt);
  
  
  
  gszspk(channelData,spks,tots,channelsLength,intervalLength);
  if (verbose==13) {
    if (pL->isz<1||pL->plen[0]!=intervalsCount) {
      printf("For verbose 13 need 1 vec of %d each\n",intervalsCount); hxe(); } 
    for (ii=0;ii<intervalsCount;ii++) pL->pv[0][ii]=(double)spks[ii];
  }
  if(incby1) for (ii=0; ii<intervalsCount; ii++) { 
    if (flag==0) {
      true_var=(spks[ii]>minspikes && diffcor[ii]>diffthresh);
    } else if (flag==1) {true_var=(spks[ii]>minspikes);
    } else if (flag==2) {true_var=(diffcor[ii]>diffthresh);
    }
    if (true_var) {
      if(seizuresCount==0 || tmpSeizure->endIntervalIndex!=ii-1)  { 
        seizuresCount++;
        AddSeizure(&tmpSeizures,tmpSeizure=AllocSeizure(ii,spks[ii]));
        AddSeizure(pSeizuresOut,AllocSeizure(ii,spks[ii]));
      } else { 
        tmpSeizure->endIntervalIndex = ii; 
        tmpSeizure->totalSpikesCount += spks[ii]; 
        tmpSeizure->endSpikesCount = spks[ii]; 
      }
    }
  } 
  else  for (ii=0; ii+1<intervalsCount; ii+=2) { 
    if (flag==0) {
      true_var=(spks[ii]>minspikes && spks[ii+1]>minspikes && diffcor[ii]>diffthresh && diffcor[ii+1]>diffthresh);
    } else if (flag==1) {true_var=(spks[ii]>minspikes && spks[ii+1]>minspikes);
    } else if (flag==2) {true_var=(diffcor[ii]>diffthresh && diffcor[ii+1]>diffthresh);
    }
    if (true_var) {
      if(seizuresCount==0 || tmpSeizure->endIntervalIndex!=ii-1)  { 
        seizuresCount++;
        AddSeizure(&tmpSeizures,tmpSeizure=AllocSeizure(ii,spks[ii]+spks[ii+1]));
        AddSeizure(pSeizuresOut,AllocSeizure(ii,spks[ii]+spks[ii+1]));
      } else { 
        tmpSeizure->endIntervalIndex = ii+1; 
        tmpSeizure->totalSpikesCount += spks[ii]+spks[ii+1]; 
        tmpSeizure->endSpikesCount = spks[ii+1]; 
      }
    }
  } 
  if(verbose==-1) {
    printf("before tmpSeizures %d:\n",tmpSeizures.count);  printSeizures(&tmpSeizures); }
  
  
  for (seizureIndex=0; seizureIndex<seizuresCount; seizureIndex++) { 
    tmpSeizure = tmpSeizures.pp[seizureIndex];
    currentSeizure = pSeizuresOut->pp[seizureIndex];
    currentSeizure->ID = seizureIndex;
    currentSeizure->totalSpikesCount = tmpSeizure->totalSpikesCount;
    if(seizureIndex < seizuresCount-1) { 
      
      
      
      stopIndex = (tmpSeizures.pp[seizureIndex+1])->startIntervalIndex * intervalLength; 
    } else stopIndex = channelsLength;
    intervalIndex = tmpSeizure->endIntervalIndex + 1; 
    
    dataIndex = intervalIndex * intervalLength; 
    if (intervalIndex >= totsLen) intervalIndex = totsLen-1;
    limit = tots[intervalIndex]/50.; 
    avgnumspikes=tmpSeizure->totalSpikesCount/\
      (Bintdur*(tmpSeizure->endIntervalIndex - tmpSeizure->startIntervalIndex+1));
    if (verbose==-2) printf("avgnumspikes: %d\n",avgnumspikes);
    upcount=upvalue=upflag=0;
    for (j=1; j<minspikes; j++) { 
      spikesCount = 0;
      for (i=0; i<samprate*Sintdur; i++, dataIndex++) { 
        if(dataIndex>=stopIndex) { 
          dataIndex = stopIndex-1; break; }
        if (dataIndex+1 >= channelsLength) break;
          value=channelData[dataIndex+1]-channelData[dataIndex];
          if(value>0) { upflag = 1; 
            upcount++; 
            upvalue+=value; 
          } else	{ 
            if (upflag && upvalue>limit && upcount>uplim) {
              if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
                sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
                if(sharp < sharpth) spikesCount++; 
              } else spikesCount++;
            }
            upflag=upcount=upvalue=0;
          }
      } 
      if (spikesCount/Sintdur < endwting*avgnumspikes) break; 
    } 
    currentSeizure->endIndex=dataIndex; 
    currentSeizure->endIntervalIndex = (int)dataIndex/intervalLength; 
  }

  
  
  
  seizuresCount = pSeizuresOut->count; 
  for (seizureIndex = 0; seizureIndex<seizuresCount; seizureIndex++) { 
    tmpSeizure = tmpSeizures.pp[seizureIndex];
    currentSeizure = pSeizuresOut->pp[seizureIndex];
    if (seizureIndex>0) { 
      
      
      
      
      stopIndex = (pSeizuresOut->pp[seizureIndex-1])->endIndex;
    }
    else stopIndex = 0;
    intervalIndex = tmpSeizure->startIntervalIndex; 
    if(intervalIndex >= totsLen) intervalIndex = totsLen-1;
    dataIndex = intervalIndex * intervalLength; 
    limit= tots[intervalIndex]/50; 
    avgnumspikes = tmpSeizure->totalSpikesCount/\
      (Bintdur*(tmpSeizure->endIntervalIndex-tmpSeizure->startIntervalIndex+1));
    upcount=0.0; upvalue=upflag=0; 
    for(j=1; j<minspikes; j++) { 
      spikesCount = 0;
      for(i=0; i<samprate*Sintdur; i++, dataIndex--) { 
	if(dataIndex<=stopIndex) { 
	  dataIndex = stopIndex + 1; break; }
        if(dataIndex+1 >= channelsLength) break;
        value=channelData[dataIndex+1]-channelData[dataIndex];
        if(value>0) { 
          upflag=1; upcount++; upvalue+=value;
        } else { 
          if (upflag && upvalue>limit && upcount>uplim) {
            if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
              sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
              if(sharp < sharpth) spikesCount++; 
            } else spikesCount++;
          }
          upflag=upcount=upvalue=0;
        }
      } 
      if (spikesCount/Sintdur<endwting*avgnumspikes) break;
    }
    currentSeizure->startIndex = dataIndex; 
    currentSeizure->startIntervalIndex = (int)dataIndex/intervalLength; 
  }
  
  for(seizureIndex = seizuresCount-1; seizureIndex>0; seizureIndex--) { 
    currentSeizure = pSeizuresOut->pp[seizureIndex-1];
    nextSeizure =    pSeizuresOut->pp[seizureIndex]; 
    if(currentSeizure->endIndex >= (nextSeizure->startIndex - mindist)) {
      
      currentSeizure->endIndex = nextSeizure->endIndex; 
      currentSeizure->totalSpikesCount += nextSeizure->totalSpikesCount;
      free(nextSeizure);     
      RemoveSeizureAt(pSeizuresOut,seizureIndex);
    }
  }
  
  if(verbose==-1) {
    printf("after tmpSeizures %d:\n",tmpSeizures.count);  printSeizures(&tmpSeizures);
    printf("after pSeizuresOut %d:\n",pSeizuresOut->count); printSeizures(pSeizuresOut); }
 STALEY_DOFREE:
  if(diffcor) free(diffcor);
  if(percentc) free(percentc);
  if(tots) free(tots);
  for(i=0; i<tmpSeizures.count; i++) { 
    if(tmpSeizures.pp[i]==0x0) printf("tmpSeizures.pp[%d]=0x0",i);
    free(tmpSeizures.pp[i]); 
  }
  FreeLSeizure(&tmpSeizures); 
  return pSeizuresOut;
}


static double gsz (double* channelData,double* diffcor,double* percentc,double* tots,\
          int channelsLength,int intervalLength, int SgroupCount,int LintCnt) {
  int intervalIndex, dataIndex, iS, intervalsCount;
  int i,j, ii, SgroupCount4, dbxi;
  double value, min_minS, max_maxS, HV, LV; 
  double diffc, totsum, avg1, sdev1, avg2, sdev2;
  double upvalue, limit, *minS, *maxS;
  double bintdur,lintdur; 
  
  intervalsCount = channelsLength / intervalLength; 
  SgroupCount4 = SgroupCount/4;  
  
  
  if((minS = (double*) calloc(SgroupCount+2,sizeof(double)))==0x0) {
    printf("getseizures ERR: out of memory!\n"); hxe();  }
  
  
  if((maxS = (double*) calloc(SgroupCount+2,sizeof(double)))==0x0){
    printf("getseizures ERR: out of memory!\n"); hxe();  }

  
  for (dbxi=0,intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++) { 
    diffcor[intervalIndex]=percentc[intervalIndex]=tots[intervalIndex]=totsum=diffc=value=0;
    
    
    for (dataIndex=intervalIndex*intervalLength, iS=0; iS<SgroupCount+2; iS++) { 
      minS[iS]= 1e22; maxS[iS]= -1e22; 
      for(i=0; i<LintCnt; i++, dataIndex++) {  
	if(dataIndex < channelsLength) value = channelData[dataIndex];
        if (minS[iS] > value)  minS[iS] = value;
        if (maxS[iS] < value)  maxS[iS] = value;
      }
    }
    if (verbose==11) {
      if (dbxi==0 && (pL->isz<2 || pL->plen[0]!=intervalsCount*(SgroupCount+2)\
          || pL->plen[1]!=intervalsCount*(SgroupCount+2))) {
        printf("For verbose 11 need 2 vecs of %d each\n",intervalsCount*(SgroupCount+2)); hxe();} 
      for (iS=0; iS<SgroupCount+2; iS++,dbxi++) {
        pL->pv[0][dbxi]=minS[iS]; pL->pv[1][dbxi]=maxS[iS];
      }
    }
    
    for (iS=0; iS<SgroupCount; iS++) {
      HV = MIN(maxS[iS], MAX(maxS[iS+1],maxS[iS+2])); 
      LV = MAX(minS[iS], MIN(minS[iS+1],minS[iS+2]));
      
      diffc += (HV - LV); 
      totsum += (maxS[iS]-minS[iS]); 
    }    
    
    diffcor[intervalIndex] = diffc;	
    percentc[intervalIndex] = diffc/totsum; 
    
    for(totsum=0, iS=0, j=0; j<4; j++)  {  
      min_minS = 1e22;  max_maxS = -1e22; 
      for(i=0; i<SgroupCount4; i++, iS++) { 
        if (max_maxS < maxS[iS])  max_maxS = maxS[iS];
        if (min_minS > minS[iS])  min_minS = minS[iS];
      }
      totsum += (max_maxS - min_minS);
    }
    tots[intervalIndex] = totsum/4;
    if (verbose==12) {
      if (dbxi==0) {
        for (ii=0;ii<3;ii++) if (pL->isz<3||pL->plen[ii]!=intervalsCount) {
          printf("For verbose 12 need 3 vecs of %d each\n",intervalsCount); hxe(); } 
        printf("Verbose 12: diffcor,percentc,tots\n");
      }
      ii=intervalIndex;
      pL->pv[0][dbxi]=diffcor[ii]; pL->pv[1][dbxi]=percentc[ii]; pL->pv[2][dbxi]=tots[ii];
      dbxi++;
    }
  } 
  
  
  
  
  
  avg1=avg2=sdev1=sdev2=0.0;
  for(intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++)  { 
    avg1 += diffcor[intervalIndex];  
    avg2 += percentc[intervalIndex]; 
    sdev1 += diffcor[intervalIndex]*diffcor[intervalIndex];   
    sdev2 += percentc[intervalIndex]*percentc[intervalIndex];
  }
  avg1 = avg1/intervalsCount; avg2 = avg2/intervalsCount; 
  sdev1 = sdev1/intervalsCount - avg1*avg1; 
  if(sdev1>0.) sdev1=sqrt(sdev1); else sdev1=avg1;
  sdev2 = sdev2/intervalsCount - avg2*avg2;
  if(sdev2>0.) sdev2=sqrt(sdev2); else sdev2=avg2;
  if (verbose>1) printf("diffcor: %g (%g,%g), percentc: %g (%g)\n",\
                        avg1,sdev1,avg1+abovebth*sdev1,avg2,sdev2);
  GSZ_DOFREE:
  if(minS) free(minS);
  if(maxS) free(maxS);
  return avg1+abovebth*sdev1;
}

#ifdef MYSPUD

int dspud (double* src, int nsrc, int lc) {
  int i, k, m, n, nqsz, nsrc, jj[UDSL], f[UDSL], lc, dsz[UDSL], nqmax, thsz, lc2, done, dbn;
  double *src, *tvec, *th, *dest[UDSL], *nq[UDNQ], *tmp, *dbx, lt, thdist;
  Object *ob, *ob2;
  void *vvd[UDSL], *vvth, *vnq[UDNQ];
  
  
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

#endif

static double gszspk (double* channelData, int* spks,\
                      double* tots, int channelsLength, int intervalLength) {
  int intervalIndex, dataIndex, intervalsCount, foundspk;
  int i,j, ii, upflag, upcount, spikesCount, uplim, dbgSpikes,didpr,cnt;
  double value, upvalue, limit, sum, sharp; 
  
  intervalsCount = channelsLength / intervalLength; 
  uplim=(int)(0.5+spkup*samprate/1e3);
  dbgSpikes=didpr=0;

  if(useavgtots) {
    sum=0.0;
    for(intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++) sum += tots[intervalIndex];
    sum /= (double) intervalsCount;
    limit= sum*spklim;
  }

  for (intervalIndex=0; intervalIndex<intervalsCount; intervalIndex++) { 
    dataIndex = intervalIndex * intervalLength;
    if(!useavgtots) limit= tots[intervalIndex]*spklim; 
    if(verbose>=15 && !didpr){ printf("limit = %g\n",limit); didpr=1; }
    upcount=upvalue=spikesCount=upflag=0;
    for (j=0; j<intervalLength; j++, dataIndex++)  {  
      if (dataIndex > channelsLength) continue; 
      value = channelData[dataIndex+1]-channelData[dataIndex];
      if(value>0) { 
        upflag=1; upcount++; upvalue+=value;
      } else { 
        foundspk=0; sharp=0.0;
        if(upflag && upvalue>limit && upcount>uplim) {
          foundspk=1;
          if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
            sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
            if(sharp > sharpth) foundspk=0; 
          }
        }
        if (foundspk) { 
          spikesCount++;  
          if(verbose>=14) { 
            if(pL->isz < 3 || pL->plen[0]<channelsLength || pL->plen[1]<channelsLength || pL->plen[2]<1) {
              printf("need at least 2 vectors of size %d for verbose==14!\n",channelsLength); hxe();
            } else { 
              if(usesharp && dataIndex-sharpoff>=0.0 && dataIndex+sharpoff<channelsLength) {
                sharp = channelData[dataIndex+(int)sharpoff]-2*channelData[dataIndex]+channelData[dataIndex-(int)sharpoff];
                if(verbose>=15) printf("spike found (x,y,upcount,upvalue,sharp)=(%d,%g,%d,%g,%g)\n",dataIndex,channelData[dataIndex],upcount,upvalue,sharp);
              } else {
                sharp=0.0;
                if(verbose>=15) printf("spike found (x,y,upcount,upvalue)=(%d,%g,%d,%g)\n",dataIndex,channelData[dataIndex],upcount,upvalue);
              }              
              pL->pv[0][dbgSpikes]=dataIndex; pL->pv[1][dbgSpikes++]=channelData[dataIndex];}
          }
        }
        upflag=upcount=upvalue=0; 
      }
    } 
    spks[intervalIndex]=spikesCount;
  }
  if(verbose>=14) pL->pv[2][0]=dbgSpikes;
  return 0.;
}



static double getseizures (void* vv) {
  int n, cnt,i; LSeizure* pSeizures;
  double *p,*totalSpikesCount,*startIndex,*endIndex;
  n = vector_instance_px(vv,&p);
  if(verbose>10) {
    if (!ifarg(4)) { printf("Use veclist for dbx with verbose>10\n");hxe();
    } else pL=AllocListVec(*hoc_objgetarg(4));
  }
  if(!(pSeizures=dgetseizures(p,n))) return 0.0;
  cnt=pSeizures->count;
  totalSpikesCount=vector_newsize(vector_arg(1),cnt);
  startIndex=vector_newsize(vector_arg(2),cnt);
  endIndex=vector_newsize(vector_arg(3),cnt);

  for(i=0;i<cnt;i++) {
    totalSpikesCount[i] =   (double)pSeizures->pp[i]->totalSpikesCount;
    startIndex[i] =         (double)pSeizures->pp[i]->startIndex;
    endIndex[i] =           (double)pSeizures->pp[i]->endIndex;
  } 
  for(i=0;i<pSeizures->count;i++) free(pSeizures->pp[i]);
  FreeLSeizure(pSeizures);
  if (pL) FreeListVec(&pL);
  return (double)cnt;
}
ENDVERBATIM

PROCEDURE install () {
  VERBATIM
  if(!installed) {
    install_vector_method("getseizures",getseizures);
  } else printf("%s\n","$Id: staley.mod,v 1.75 2010/04/30 18:59:54 samn Exp $");
  ENDVERBATIM
  installed=1
}