NEURON {
  SUFFIX intfsw
  GLOBAL INSTALLED
  GLOBAL verbose
  GLOBAL edgefuncid 
}

PARAMETER {
  INSTALLED=0
  verbose=0
  edgefuncid=0
}

VERBATIM
#include "misc.h"

typedef struct {
  int isz;
  int imaxsz;
  double* p;  
} myvec;

myvec* allocmyvec (int maxsz){
  myvec* pv = (myvec*)malloc(sizeof(myvec));
  if(!pv) return 0x0;
  pv->isz=0;
  pv->imaxsz=maxsz;
  pv->p=(double*)malloc(sizeof(double)*maxsz);
  if(!pv->p) { free(pv); return 0x0; }
  return pv;
}

int freemyvec (myvec** pps) {
  if(!pps || !pps[0]) return 0;
  myvec* ps = pps[0];
  if(ps->p)free(ps->p);
  free(ps);
  pps[0]=0x0;
  return 1;
}

double popmyvec (myvec* pv) {
  if(pv->isz<1) {
    printf("popmyvec ERRA: can't pop empty stack!\n");
    return 0.0;
  }
  double d = pv->p[pv->isz-1]; pv->isz--;
  return d;
}

void popallmyvec (myvec* pv) {
  pv->isz=0;
}

double pushmyvec (myvec* ps,double d) {
  if(ps->isz==ps->imaxsz) {
    printf("pushmyvec realloc\n");
    ps->imaxsz*=2;
    ps->p=(double*)realloc(ps->p,sizeof(double)*ps->imaxsz);
    if(!ps->p){ printf("pushmyvec ERRA: myvec out of memory %d!!\n",ps->imaxsz); return 0.0; }
  }
  ps->p[ps->isz++]=d; 
  return 1.0;  
}

double appendmyvec (myvec* ps,double d) {
  return pushmyvec(ps,d);
}

typedef struct myqnode_ {
  struct myqnode_* pnext;  
  struct myqnode_* pprev;
  int dd;
} myqnode;

myqnode* allocmyqnode() {
  myqnode* p = (myqnode*)malloc(sizeof(myqnode));
  p->pnext=0x0;
  p->pprev=0x0;
  return p;
}

typedef struct {
  myqnode* pfront;
  myqnode* pback;
} myq;

myq* allocmyq() {
  myq* pq = (myq*)malloc(sizeof(myq));
  pq->pfront = pq->pback = 0x0;
  return pq;
}

int freemyq(myq** ppq) {
  myq* pq = *ppq;
  myqnode* ptmp=pq->pback;
  while(pq->pback){
    if(pq->pback->pprev==0x0){
      free(pq->pback);
      pq->pback=0x0;
      pq->pfront=0x0;
      break;
    } else {
      ptmp=pq->pback->pprev;
      free(pq->pback);    
    }
  }
  free(pq);
  ppq[0]=0;
  return 1;
}

int printfrontmyq (myq* pq) {
  if(pq && pq->pfront) {
    printf("front=%d  ",pq->pfront->dd);
    return 1;
  }
  printf("printfrontmyq ERRA: empty front!\n");
  return 0;
}

int printbackmyq (myq* pq) {
  if(pq && pq->pback) {
    printf("back=%d  ",pq->pback->dd);
    return 1;
  }
  printf("printbackmyq ERRA: empty back!\n");
  return 0;
}

int printmyq (myq* pq, int backwards) {
  if(pq){
    int i=0;
    if(backwards){
      myqnode* pnode = pq->pback;
      while(pnode){
        printf("val %d from back = %d\n",i++,pnode->dd);
        pnode = pnode->pprev;
      }
    } else {
      myqnode* pnode = pq->pfront;
      while(pnode){
        printf("val %d from front = %d\n",i++,pnode->dd);
        pnode = pnode->pnext;
      }
    }
    return 1;
  }
  printf("printmyq ERRA: null pointer!\n");
  return 0;
}

int enqmyq (myq* pq,int d) {
  if(pq->pfront==pq->pback) {
    if(!pq->pfront){
      pq->pfront = allocmyqnode();
      pq->pback = pq->pfront;
      pq->pfront->dd=d;
    } else {
      pq->pback = allocmyqnode();
      pq->pback->dd=d;
      pq->pback->pprev = pq->pfront;
      pq->pfront->pnext = pq->pback;
    }
  } else {
    myqnode* pnew = allocmyqnode();
    pnew->dd = d;
    pq->pback->pnext = pnew; 
    pnew->pprev = pq->pback;
    pq->pback = pnew;
  }
  return 1;
}

int emptymyq (myq* pq) {
  if(pq->pfront==0x0) return 1;
  return 0;
}

int deqmyq (myq* pq) {
  if(pq->pfront == pq->pback){
    if(!pq->pfront){
      printf("deqmyq ERRA: can't deq empty q!\n");
      return -1.0;
    } else {
      int d = pq->pfront->dd;
      free(pq->pfront);
      pq->pfront=pq->pback=0x0;
      return d;
    }
  } else {
    myqnode* tmp = pq->pfront;
    int d = tmp->dd;
    pq->pfront = pq->pfront->pnext;
    pq->pfront->pprev = 0x0;
    free(tmp);
    return d;
  }
}

ENDVERBATIM

FUNCTION testmystack () {
VERBATIM
  myvec* pv = allocmyvec(10);
  printf("created stack with sz %d\n",pv->imaxsz);
  int i;
  for(i=0;i<pv->imaxsz;i++) {
    double d = 41.0 * (i%32) + rand()%100;
    printf("pushing %g onto stack of sz %d\n",d,pv->isz);
    pushmyvec(pv,d);
  }
  printf("test stack realloc by pushing 123.0\n");
  pushmyvec(pv,123.0);
  printf("stack now has %d elements, %d maxsz. contents:\n",pv->isz,pv->imaxsz);
  for(i=0;i<pv->isz;i++)printf("s[%d]=%g\n",i,pv->p[i]);
  printf("popping %d elements. contents:\n",pv->isz);
  while(pv->isz){
    double d = popmyvec(pv);
    printf("popped %g, new sz = %d\n",d,pv->isz);
  }
  printf("can't pop stack now, empty test: ");
  popmyvec(pv);
  freemyvec(&pv);
  printf("freed stack\n");
  return 1.0;
ENDVERBATIM
}

FUNCTION testmyq () {
VERBATIM
  myq* pq = allocmyq();
  printf("created q, empty = %d\n",emptymyq(pq));
  printf("enqueing 10 values:\n");
  int i;
  for(i=0;i<10;i++){
    int d = 41 * (i%32) + rand()%252;
    printf("enqueuing %d...",d);
    enqmyq(pq,d);
    printfrontmyq(pq);
    printbackmyq(pq); printf("\n");
  }
  printf("printing q in forwards order:\n");
  printmyq(pq,0);
  printf("printing q in backwards order:\n");
  printmyq(pq,1);
  printf("testing deq:\n");
  while(!emptymyq(pq)){
    printf("b4 deq: ");
    printfrontmyq(pq); 
    printbackmyq(pq); printf("\n");
    int d = deqmyq(pq);
    printf("dequeued %d\n",d);
    printf("after deq: ");
    printfrontmyq(pq); 
    printbackmyq(pq); printf("\n");
  }
  freemyq(&pq);
  printf("freed myq\n");
  return 1.0;
ENDVERBATIM
}


VERBATIM


static double copynz (void* vv) {
  double* pV;
  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;
  int iStartIDx = 0, iEndIDx = n - 1;
  if(ifarg(2)){
    iStartIDx = (int)*getarg(1);
    iEndIDx = (int) *getarg(2);
  }
  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n
                         || iStartIDx<0    || iEndIDx < 0){
    printf("copynz ERRA: invalid indices start=%d end=%d size=%d\n",iStartIDx,iEndIDx,n);
    return -1.0;
  }

  double* pVal,*pBin;

  if(vector_arg_px(1,&pVal)!=n || vector_arg_px(2,&pBin)!=n){
    printf("copynz ERRB: vec args must have size %d!",n);
    return 0.0;
  }

  int iOutSz = 0;
  for(idx=iStartIDx;idx<=iEndIDx;idx++){
    if(pBin[idx]){
      pV[iOutSz++]=pVal[idx];
    }
  }

  vector_resize(pV,iOutSz);

  return (double)iOutSz;
}


static double nnmeandbl (double* p,int iStartIDX,int iEndIDX) {
  int iCount=0,idx=0;
  double dSum = 0.0;
  for(idx=iStartIDX;idx<=iEndIDX;idx++){
    if(p[idx]>=0.0){
      dSum+=p[idx];
      iCount++;
    }
  }
  if(iCount>0) return dSum / iCount;
  return -1.0;
} 


static double gzmeandbl (double* p,int iStartIDX,int iEndIDX) {
  int iCount=0,idx=0;
  double dSum = 0.0;
  for(idx=iStartIDX;idx<=iEndIDX;idx++){
    if(p[idx]>0.0){
      dSum+=p[idx];
      iCount++;
    }
  }
  if(iCount>0) return dSum / iCount;
  return -1.0;
}


static double gzmean (void* vv) {
  double* pV;
  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;
  int iStartIDx = 0, iEndIDx = n - 1;
  if(ifarg(2)){
    iStartIDx = (int)*getarg(1);
    iEndIDx = (int) *getarg(2);
  }
  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n
                         || iStartIDx<0    || iEndIDx < 0){
    printf("gzmean ERRA: invalid indices start=%d end=%d size=%d\n",iStartIDx,iEndIDx,n);
    return -1.0;
  }
  return gzmeandbl(pV,iStartIDx,iEndIDx);
}



static double nnmean (void* vv) {
  double* pV;
  int n = vector_instance_px(vv,&pV) , iCount = 0 , idx=0;
  int iStartIDx = 0, iEndIDx = n - 1;
  if(ifarg(2)){
    iStartIDx = (int)*getarg(1);
    iEndIDx = (int) *getarg(2);
  }
  if(iEndIDx < iStartIDx || iStartIDx >= n || iEndIDx >= n
                         || iStartIDx<0    || iEndIDx < 0){
    printf("nnmean ERRA: invalid indices start=%d end=%d size=%d\n",iStartIDx,iEndIDx,n);
    return -1.0;
  }
  return nnmeandbl(pV,iStartIDx,iEndIDx);
}
ENDVERBATIM


FUNCTION GetCCR () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCC ERRA: problem initializing first arg!\n");
    return 0.0;
  }

  int iCells = pList->isz;
  if(iCells<2){
    printf("GetCC ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  int* pNeighbors = (int*)calloc(iCells,sizeof(int));
  int i = 0, iNeighbors = 0;
  if(!pNeighbors){
    printf("GetCCR ERRE: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }  

  
  double* pCC; 
  int iVecSz = vector_arg_px(2,&pCC);
  if(!pCC || iVecSz < iCells){
    printf("GetCCR ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pCC,0,sizeof(double)*iVecSz);

  
  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,
      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetCCR ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       free(pNeighbors);
       return 0.0;
  }

  double dSubsamp = ifarg(5)?*getarg(5):1.0;
  if(dSubsamp<0.0 || dSubsamp>1.0){
    printf("GetCCR ERRH: invalid subsamp = %g , must be btwn 0 and 1\n",dSubsamp);
    FreeListVec(&pList);
    free(pNeighbors);
    return 0.0;
  }

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  double* pUse = 0; 
  
  if(dSubsamp<1.0){ 
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  
  int myID;

  int* pNeighborID = (int*)calloc(iCells,sizeof(int));

  if( verbose > 0 ) printf("searching from id: ");

  for(myID=0;myID<iCells;myID++) pCC[myID]=-1.0; 

  for(myID=iStartID;myID<=iEndID;myID++){

    if(verbose > 0 && myID%1000==0)printf("%d ",myID);

    
    if(pUse && pUse[myID]>=dSubsamp) continue;

    int idx = 0, youID = 0, youKidID=0 , iNeighbors = 0;

    
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID>=iStartID && youID<=iEndID){
        pNeighbors[youID]=1;      
        pNeighborID[iNeighbors++]=youID;
      }
    }

    if(iNeighbors < 2){
      for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
      continue;
    }

    int iConns = 0 ; 
  
    
    for(i=0;i<iNeighbors;i++){
      if(!pNeighbors[pNeighborID[i]])continue;
      youID=pNeighborID[i];
      for(idx=0;idx<pLen[youID];idx++){
        youKidID=pLV[youID][idx];
        if(youKidID >= iStartID && youKidID <= iEndID && pNeighbors[youKidID]){
          iConns++;
        }
      }
    }
    pCC[myID]=(double)iConns/((double)iNeighbors*(iNeighbors-1));
    for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
  }
 
  free(pNeighborID);
  free(pNeighbors);
  FreeListVec(&pList);
  if(pUse)free(pUse);

  if( verbose > 0 ) printf("\n");

  return  1.0;
  ENDVERBATIM
}




FUNCTION GetCentrality () {
  VERBATIM

  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCentrality ERRA: problem initializing first arg!\n");
    return 0.0;
  }

  int iCells = pList->isz;
  if(iCells<2){
    printf("GetCentrality ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  double* pCE; 
  int iVecSz = vector_arg_px(2,&pCE);
  if(!pCE || iVecSz < iCells){
    printf("GetCCR ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pCE,0,sizeof(double)*iVecSz);

  double dSubsamp = ifarg(3)?*getarg(3):1.0;
  if(dSubsamp<0.0 || dSubsamp>1.0){
    printf("GetCCR ERRH: invalid subsamp = %g , must be btwn 0 and 1\n",dSubsamp);
    FreeListVec(&pList);
    return 0.0;
  }

  unsigned int iSeed = ifarg(4)?(unsigned int)*getarg(4):INT_MAX-109754;

  double* pUse = 0; 
  
  if(dSubsamp<1.0){ 
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  int s,w,T,v,idx;

  myvec* S = allocmyvec(iCells*2);
  myvec** P = (myvec**)malloc(sizeof(myvec*)*iCells);
  myvec* d = allocmyvec(iCells);
  myvec* sigma = allocmyvec(iCells);
  myvec* di = allocmyvec(iCells);
  for(w=0;w<iCells;w++) P[w]=allocmyvec(iCells);
  for(s=0;s<iCells;s++){
    if(verbose && s%100==0) printf("s=%d\n",s);
    S->isz=0;
    for(w=0;w<iCells;w++) P[w]->isz=0;
    for(T=0;T<iCells;T++) sigma->p[T]=0; sigma->p[s]=1;
    for(T=0;T<iCells;T++) d->p[T]=-1; d->p[s]=0;
    myq* Q = allocmyq();
    enqmyq(Q,s);
    while(!emptymyq(Q)){
      v = deqmyq(Q);
      pushmyvec(S,v);
      for(idx=0;idx<pLen[v];idx++){
        w = (int) pLV[v][idx];
        if(d->p[w]<0){
          enqmyq(Q,w);
          d->p[w] = d->p[v] + 1;
        }
        if(d->p[w] == d->p[v] + 1){
          sigma->p[w] = sigma->p[w] + sigma->p[v];
          appendmyvec(P[w],v);
        }
      }
    }
    freemyq(&Q);
    for(v=0;v<iCells;v++) di->p[v]=0;
    while(S->isz){
      w = popmyvec(S);
      for(idx=0;idx<P[w]->isz;idx++){
        v=P[w]->p[idx];
        di->p[v] = di->p[v] + (sigma->p[v]/sigma->p[w])*(1.0+di->p[w]);
      }
      if(w!=s) pCE[w] = pCE[w] + di->p[w];
    }
  }

  int N = 0;
  for(s=0;s<iCells;s++) if(pLen[s]) N++;
  if(N>2){
    double scale = 1.0/( (N-1.0)*(N-2.0) );
    for(v=0;v<iCells;v++) if(pLen[v]) pCE[v] *= scale;
  }
  
CEFREE:
  freemyvec(&S);
  for(w=0;w<iCells;w++) freemyvec(&P[w]);
  free(P);
  freemyvec(&d);
  freemyvec(&sigma);
  freemyvec(&di);
  if(pUse)free(pUse);  
  return 1.0;

  ENDVERBATIM
}






FUNCTION GetCC () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCC ERRA: problem initializing first arg!\n");
    return -1.0;
  }

  int iCells = pList->isz;
  if(iCells<2){
    printf("GetCC ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return -1.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  int* pNeighbors = (int*)calloc(iCells,sizeof(int));
  int i = 0, iNeighbors = 0;
  if(!pNeighbors){
    printf("GetCC ERRE: out of memory!\n");
    FreeListVec(&pList);
    return -1.0;
  }  

  
  int myID = (int) *getarg(2);
  if(myID < 0 || myID >= iCells){
    printf("GetCC ERRF: invalid id = %d\n",myID);
    FreeListVec(&pList);
    free(pNeighbors);
    return -1.0;
  }

  
  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,
      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetCC ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       free(pNeighbors);
       return -1.0;
     }

  int idx = 0, iDist = 1 , youID = 0, youKidID=0;

  int* pNeighborID = (int*)calloc(iCells,sizeof(int));

  
  for(idx=0;idx<pLen[myID];idx++){
    youID = pLV[myID][idx];
    if(youID>=iStartID && youID<=iEndID){
      pNeighbors[youID]=1;      
      pNeighborID[iNeighbors++]=youID;
    }
  }

  if(iNeighbors < 2){
    FreeListVec(&pList);
    free(pNeighbors);
    return -1.0;
  }

  int iConns = 0; 

  
  for(i=0;i<iNeighbors;i++){
    if(!pNeighbors[pNeighborID[i]])continue;
    youID=pNeighborID[i];
    for(idx=0;idx<pLen[youID];idx++){
      youKidID=pLV[youID][idx];
      if(youKidID >= iStartID && youKidID <= iEndID && pNeighbors[youKidID]){
        iConns++;
      }
    }
  }
 
  free(pNeighborID);
  free(pNeighbors);
  FreeListVec(&pList);

  return  (double)iConns/((double)iNeighbors*(iNeighbors-1));
  
  ENDVERBATIM
}








FUNCTION CountNeighborsR () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("CountNeighborsR ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("CountNeighborsR ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("CountNeighborsR ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);

  
  int myID = (int) *getarg(3);
  if(myID < 0 || myID >= iCells){
    printf("CountNeighborsR ERRF: invalid id = %d\n",myID);
    FreeListVec(&pList);
    return 0.0;
  }

  
  int iStartID = (int)*getarg(3),
      iEndID =   (int)*getarg(4),
      iSearchDegree =    (int)*getarg(5);

  double dSubsamp = ifarg(6)?*getarg(6):1.0;

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("CountNeighborsR ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       return 0.0;
     }

  
  if(iSearchDegree<=0){
    printf("CountNeighborsR ERRI: invalid searchdegree=%d\n",iSearchDegree);
    FreeListVec(&pList);
    return 0.0;
  }

  
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("CountNeighborsR ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0, iMatches = 0;

  double* pVDTmp = 0, dgzt = 0.0; 
  int* pTmp = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ 
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  if( verbose > 0 ) printf("searching from id: ");

  pVDTmp = (double*)calloc(iCells,sizeof(double));
  pTmp = (int*)calloc(iCells,sizeof(int)); 

  for(myID=iStartID;myID<=iEndID;myID++){

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    
    if(pUse && pUse[myID]>=dSubsamp) continue;

    iMatches = 0;

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;

    
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID>=iStartID && youID<=iEndID && !pVDTmp[youID]){
        pVDTmp[youID]=(double)iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    if(iSearchDegree == iDist){
      pVD[myID] = iCheckSz;
      for(idx=0;idx<iCheckSz;idx++) pVDTmp[pCheck[idx]]=0; 
      continue;
    }

    pVDTmp[myID]=1;

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    
    while(iCheckSz>0 && iDist<=iSearchDegree){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(youKidID >= iStartID && youKidID <=iEndID && !pVDTmp[youKidID]){ 
            pTmp[iTmpSz++] = youKidID; 
            pVDTmp[youKidID]=(double)iDist; 
          }
        }
      }
      iCheckSz = iTmpSz;
      
      if(iSearchDegree == iDist){
        pVD[myID] = iCheckSz;
        memset(pVDTmp,0,sizeof(double)*iCells); 
        break;
      } 

      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }
  }

  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp); free(pTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
  ENDVERBATIM
}


VERBATIM
double maxval(double* p,int sz)
{
  double dmax = p[0];
  int i = 1;
  for(;i<sz;i++) if(p[i]>dmax) dmax = p[i];
  return dmax;
}

double weightdelaydist(double w,double d)
{
  if(w < 0)
    return -w/d;
  if(w > 0)
    return d/w;
  return DBL_MAX; 
}

double weightdist(double w,double d)
{
  if(w < 0)
    return -w;
  if(w > 0)
    return 1/w;
  return DBL_MAX; 
}

double delaydist(double w,double d)
{
  return d;
}

void printedgefunc(int id)
{
  switch(id){
    case 0:
     printf("weightdelaydist\n");
     break;
    case 1:
     printf("weightdist\n");
     break;
    case 2:
     printf("delaydist\n");
     break;
    default:
     printf("unknown!\n");
     break;
  }
}

ENDVERBATIM


FUNCTION predgefunc () {
  VERBATIM
  int i;
  if(ifarg(1)){ printf("%d=",(int)*getarg(1)); printedgefunc((int)*getarg(1)); printf("\n"); }    
  else for(i=0;i<3;i++){ printf("%d=",i); printedgefunc(i); printf("\n"); }
  return 0.0;
  ENDVERBATIM
}








FUNCTION GetWPath () {
  VERBATIM

  double* ppre = 0, *ppo = 0, *pwght = 0, *pdel = 0, *pout = 0;
  int iSz,iTmp,i,j,k,l;
  void* voi;

  iSz = vector_arg_px(1,&ppre);

  if(iSz < 1)
  { printf("GetWPath ERRO: invalid size for presynaptic ID Vector (arg 1) %d!\n",iSz);
    return -666.666;
  }

  if( (iTmp=vector_arg_px(2,&ppo)) != iSz)
  { printf("GetWPath ERRA: incorrectly sized postsynaptic ID Vector (arg 2) %d %d!",iSz,iTmp);
    return -666.666;
  }
  if( (iTmp=vector_arg_px(3,&pwght)) != iSz)
  { printf("GetWPath ERRB: incorrectly sized weight Vector (arg 3) %d %d!\n",iSz,iTmp);
    return -666.666;
  }
  if( (iTmp=vector_arg_px(4,&pdel)) != iSz)
  { printf("GetWPath ERRC: incorrectly sized delay Vector (arg 4) %d %d!\n",iSz,iTmp);
    return -666.666;
  }

  int maxid = maxval(ppre,iSz);

  iTmp = maxval(ppo,iSz);
  if(iTmp > maxid) maxid=iTmp;

  voi = vector_arg(5);

  if( (iTmp=vector_arg_px(5,&pout))!= maxid+1 && 0)
  { printf("GetWPath ERRD: incorrectly sized output Vector (arg 5) %d %d!\n",maxid+1,iTmp);
    return -666.666;
  }
  memset(pout,0,sizeof(double)*iTmp);

  double (*EdgeFunc)(double,double) = &weightdelaydist;
  int iEdgeFuncID = (int)edgefuncid; 
  if(iEdgeFuncID < 0 || iEdgeFuncID > 2)
  {  printf("GetWPath ERRK: invalid edgedfunc id %d!\n",iEdgeFuncID);
     return -666.666;
  } else if(iEdgeFuncID == 1) EdgeFunc = &weightdist;
    else if(iEdgeFuncID == 2) EdgeFunc = &delaydist;
  if(verbose) printedgefunc(iEdgeFuncID);

 int** adj = (int**) calloc(maxid+1,sizeof(int*));
 if(!adj)
 { printf("GetWPath ERRE: out of memory!\n");
   return -666.666;
 }

 
 
 
 double** pdist = (double**) calloc(maxid+1,sizeof(double*));

 int* pcounts = (int*) calloc(maxid+1,sizeof(int));

 
 for(i=0;i<iSz;i++)
 { 
   if(i+1<iSz && ppre[i]==ppre[i+1] && ppo[i]==ppo[i+1])
   { if(verbose>1) printf("first check double synapse i=%d\n",i);
     while(1)
     { if(i+1>=iSz) break;
       if(ppre[i]!=ppre[i+1] || ppo[i]!=ppo[i+1])
       { 
         i--;
         break;
       }
       i++; 
     }      
   }
   pcounts[(int)ppre[i]]++;    
 }

 
 for(i=0;i<maxid+1;i++){
   if(pcounts[i]){
     adj[i] = (int*)calloc(pcounts[i],sizeof(int));
     pdist[i] = (double*)calloc(pcounts[i],sizeof(double));
   }
 }

 
 int* pidx = (int*) calloc(maxid+1,sizeof(int));

 
 for(i=0;i<iSz;i++)
 { int myID = (int)ppre[i];
   if(!pcounts[myID]) continue;
   double dist = EdgeFunc(pwght[i],pdel[i]);
   j=i; 
   
   if(i+1<iSz && ppre[i]==ppre[i+1] && ppo[i]==ppo[i+1])
   { if(verbose>1) printf("check double syn i=%d\n",i);
     while(1)
     { if(i+1>=iSz) break;
       if(ppre[i]!=ppre[i+1] || ppo[i]!=ppo[i+1])
       { 
         i--;
         break;
       }
       if(j!=i) 
         dist += EdgeFunc(pwght[i],pdel[i]);
       i++; 
     }      
   }
   pdist[myID][pidx[myID]] = dist;
   adj[myID][pidx[myID]] = ppo[i];
   pidx[myID]++;
 }

 free(pidx);

 
 
 double* d = (double*) malloc( (maxid+1)*sizeof(double) ); 
 for(i=0;i<=maxid;i++)
 { if(i%100==0) printf("%d ",i);
   if(!pcounts[i])continue;
   for(j=0;j<=maxid;j++) d[j] = DBL_MAX; 
   d[i] = 0.0; 
   int changed = 0;
   for(j=0;j<maxid;j++)
   { changed=0;
     for(k=0;k<=maxid;k++) 
     { for(l=0;l<pcounts[k];l++) 
       {  if(d[adj[k][l]] > d[k] + pdist[k][l]){
            d[adj[k][l]] = d[k] + pdist[k][l];
            changed=1;
          }
       }
     }
     if(!changed){ if(verbose>1) printf("early term @ j=%d\n",j); break; }
   }












   double avg = 0.0;   
   int N = 0;
   for(j=0;j<=maxid;j++)
   { if(j!=i && d[j] < DBL_MAX)
     { avg += d[j];
       N++;
     }
   }
   if(N) pout[i] = avg / (double) N;
 }

 free(d);

 
 free(pcounts);

 for(i=0;i<=maxid;i++){
   if(adj[i]) free(adj[i]);
   if(pdist[i]) free(pdist[i]);
 }

 free(adj);
 free(pdist);

 vector_resize(voi,maxid+1); 

 return gzmeandbl(pout,0,maxid);

 ENDVERBATIM
}








FUNCTION GetPathR () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPathEV ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPathEV ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetPathEV ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);

  
  int iStartID = ifarg(3) ? (int)*getarg(3) : 0,
      iEndID = ifarg(4) ? (int)*getarg(4) : iCells - 1,
      iMaxDist = ifarg(5)? (int)*getarg(5): -1;

  double dSubsamp = ifarg(6)?*getarg(6):1.0;

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetPathEV ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       return 0.0;
     }

  
  if(iMaxDist==0){
    printf("GetPathEV ERRI: invalid maxdist=%d\n",iMaxDist);
    FreeListVec(&pList);
    return 0.0;
  }

  
  int* pCheck;
  pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPathEV ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  double* pVDTmp = 0, dgzt = 0.0; 
  int* pTmp = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ 
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching from id: ");

  pVDTmp = (double*)calloc(iCells,sizeof(double));

  int myID;

  for(myID=iStartID;myID<=iEndID;myID++){

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    
    if(pUse && pUse[myID]>=dSubsamp) continue;

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;

    pVDTmp[myID]=1;

    
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID>=iStartID && youID<=iEndID && !pVDTmp[youID]){
        pVDTmp[youID]=(double)iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    
    while(iCheckSz>0 && (iMaxDist==-1 || iDist<=iMaxDist)){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(youKidID >= iStartID && youKidID <=iEndID && !pVDTmp[youKidID]){ 
            pTmp[iTmpSz++] = youKidID; 
            pVDTmp[youKidID]=(double)iDist;
          }
        }
      }
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }

    pVDTmp[myID]=0.0; 
    if((dgzt=gzmeandbl(pVDTmp,iStartID,iEndID))>0.0) pVD[myID]=dgzt;

    memset(pVDTmp,0,sizeof(double)*iCells);
  }
  
  free(pTmp);
  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
  ENDVERBATIM
}








FUNCTION GetCCSubPop () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetCCSubPop ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetCCSubPop ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  int* pNeighbors = (int*)calloc(iCells,sizeof(int));
  int i = 0, iNeighbors = 0;
  if(!pNeighbors){
    printf("GetCCSubPop ERRE: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }  

  
  double* pCC; 
  int iVecSz = vector_arg_px(2,&pCC);
  if(!pCC || iVecSz < iCells){
    printf("GetCCSubPop ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pCC,0,sizeof(double)*iVecSz);

  double* pStart,  
          *pEnd;   

  if( vector_arg_px(3,&pStart) < iCells || vector_arg_px(4,&pEnd) < iCells){
    printf("GetCCSubPop ERRF: arg 3,4 must be Vectors with size >= %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }
  double dSubsamp = ifarg(5)?*getarg(5):1.0;

  unsigned int iSeed = ifarg(6)?(unsigned int)*getarg(6):INT_MAX-109754;

  double* pUse = 0; 
  
  if(dSubsamp<1.0){ 
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  
  int myID;

  int* pNeighborID = (int*)calloc(iCells,sizeof(int));

  if( verbose > 0 ) printf("searching from id: ");

  for(myID=0;myID<iCells;myID++) pCC[myID]=-1.0; 

  for(myID=0;myID<iCells;myID++){

    if(!pStart[myID]) continue;

    if(verbose > 0 && myID%1000==0)printf("%d ",myID);

    
    if(pUse && pUse[myID]>=dSubsamp) continue;

    int idx = 0, youID = 0, youKidID=0 , iNeighbors = 0;

    
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(pEnd[youID] && !pNeighbors[youID]){
        pNeighbors[youID]=1;      
        pNeighborID[iNeighbors++]=youID;
      }
    }

    if(iNeighbors < 2){
      for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
      continue;
    }

    int iConns = 0 ; 
  
    
    for(i=0;i<iNeighbors;i++){
      if(!pNeighbors[pNeighborID[i]])continue;
      youID=pNeighborID[i];
      for(idx=0;idx<pLen[youID];idx++){
        youKidID=pLV[youID][idx];
        if(pEnd[youKidID] && pNeighbors[youKidID]){
          iConns++;
        }
      }
    }
    pCC[myID]=(double)iConns/((double)iNeighbors*(iNeighbors-1));
    for(i=0;i<iNeighbors;i++)pNeighbors[pNeighborID[i]]=0;
  }
 
  free(pNeighborID);
  free(pNeighbors);
  FreeListVec(&pList);
  if(pUse)free(pUse);

  if( verbose > 0 ) printf("\n");

  return  1.0;

  ENDVERBATIM
}




FUNCTION GetRecurCount () {
  VERBATIM
  ListVec* pList;
  int iCells,*pLen,iFromSz,iThruSz,idx,myID,youID,jdx,iCheckSz,*pVisited,*pCheck;
  double **pLV,*pFrom,*pThru,*pR;

  pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetRecurCount ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  iCells = pList->isz; 
  if(iCells < 2){
    printf("GetRecurCount ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  pLV = pList->pv;
  pLen = pList->plen;

  pFrom=pThru=0;
  iFromSz = vector_arg_px(3,&pFrom); iThruSz = vector_arg_px(4,&pThru);
  
  if( iFromSz <= 0 || iThruSz <= 0){
    printf("GetRecurCount ERRF: arg 3,4 bad (fromsz,thrusz)=(%d,%d)\n",iFromSz,iThruSz);
    FreeListVec(&pList);
    return 0.0;
  }

  pVisited = (int*)calloc(iCells,sizeof(int));

  pCheck = (int*)malloc(sizeof(int)*iCells);

  pR = vector_newsize(vector_arg(2),iCells);
  memset(pR,0,sizeof(double)*iCells); 

  for(myID=0;myID<iCells;myID++) {
    if(!pFrom[myID]) continue;
    iCheckSz = 0; 
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(!pThru[youID] || pVisited[youID]) continue;
      pCheck[iCheckSz++]=youID;
      pVisited[youID]=1;
    }
    for(idx=0;idx<iCheckSz;idx++) {
      youID = pCheck[idx];
      for(jdx=0;jdx<pLen[youID];jdx++) {
        if(pLV[youID][jdx]==myID) pR[myID]++;
      }
    }
    memset(pVisited,0,sizeof(int)*iCells);
  }
  

  free(pCheck);
  FreeListVec(&pList);  
  free(pVisited);

  if( verbose > 0) printf("\n");

  return 1.0;

  ENDVERBATIM
}







FUNCTION GetPairDist () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPairDist ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPairDist ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  double* pFrom = 0, *pTo = 0;
  int iFromSz = vector_arg_px(3,&pFrom) , iToSz = vector_arg_px(4,&pTo);
  
  if( iFromSz <= 0 || iToSz <= 0){
    printf("GetPairDist ERRF: arg 3,4 bad (fromsz,tosz)=(%d,%d)\n",iFromSz,iToSz);
    FreeListVec(&pList);
    return 0.0;
  }

  int iMinSz = iFromSz * iToSz;

  
  double* pVD; 
  pVD = vector_newsize(vector_arg(2),iMinSz);
  memset(pVD,0,sizeof(double)*iMinSz); 

  
  int* pCheck;
  pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPairDist ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  int* pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching from id: ");

  int myID , iOff = 0 , kdx = 0;

  int* pVisited = (int*)calloc(iCells,sizeof(int)); 
  int* pUse = (int*)calloc(iCells,sizeof(int)); 
  int* pMap = (int*)calloc(iCells,sizeof(int)); 
  for(idx=0;idx<iToSz;idx++){
    pUse[(int)pTo[idx]]=1;
    pMap[(int)pTo[idx]]=idx;
  }

  for(kdx=0;kdx<iFromSz;kdx++,iOff+=iToSz){
    myID=pFrom[kdx];
    if(verbose > 0 && myID%100==0)printf("%d\n",myID);

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;
      
    
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(pUse[youID]) pVD[ iOff + pMap[youID]  ] = 1; 
      if(!pVisited[youID]){ 
        pCheck[iCheckSz++]=youID;
        pVisited[youID]=1;
      }
    }

    iTmpSz = 0;  jdx=0;
      
    iDist++;
  
    
    while(iCheckSz>0){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(pUse[youKidID] && !pVD[iOff + pMap[youKidID]])
            pVD[iOff + pMap[youKidID]] = iDist; 
          if(!pVisited[youKidID]){ 
            pTmp[iTmpSz++] = youKidID; 
            pVisited[youKidID]=1;
          }
        }
      }
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }
    memset(pVisited,0,sizeof(int)*iCells);
  }
  
  free(pTmp);
  free(pCheck);
  FreeListVec(&pList);  
  free(pUse);
  free(pMap);
  free(pVisited);

  if( verbose > 0) printf("\n");

  return 1.0;
  ENDVERBATIM
}










FUNCTION GetPathSubPop () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPathEV ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPathEV ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetPathEV ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);

  double* pStart,  
          *pEnd;   

  if( vector_arg_px(3,&pStart) < iCells || vector_arg_px(4,&pEnd) < iCells){
    printf("GetPathSubPop ERRF: arg 3,4 must be Vectors with size >= %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }
  double dSubsamp = ifarg(5)?*getarg(5):1.0;

  int bSelfLoop = ifarg(6)?(int)*getarg(6):0;

  unsigned int iSeed = ifarg(7)?(unsigned int)*getarg(7):INT_MAX-109754;

  
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPathEV ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  double  dgzt = 0.0; 
  int* pTmp = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ 
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching from id: ");

  int* pVDTmp = (int*)calloc(iCells,sizeof(int)) , myID;

  for(myID=0;myID<iCells;myID++){

    if(!pStart[myID]) continue;

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    
    if(pUse && pUse[myID]>=dSubsamp) continue;

    unsigned long int iSelfLoopDist = LONG_MAX;
    int bFindThisSelfLoop = bSelfLoop && pEnd[myID]; 

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0;

    pVDTmp[myID]=1;

    
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(bFindThisSelfLoop && youID==myID && iDist<iSelfLoopDist) iSelfLoopDist = iDist; 
      if(!pVDTmp[youID]){
        pVDTmp[youID]=iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    
    while(iCheckSz>0){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(bFindThisSelfLoop && youKidID==myID && iDist<iSelfLoopDist) iSelfLoopDist = iDist; 
          if(!pVDTmp[youKidID]){ 
            pTmp[iTmpSz++] = youKidID; 
            pVDTmp[youKidID]=iDist;
          }
        }
      }
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }

    if(bFindThisSelfLoop && iSelfLoopDist<LONG_MAX){
      pVDTmp[myID] = iSelfLoopDist;
    } else {
      pVDTmp[myID]=0; 
    }
    pVD[myID] = 0.0;
    int N = 0; 
    for(idx=0;idx<iCells;idx++){
      if(pEnd[idx] && pVDTmp[idx]){
        pVD[myID] += pVDTmp[idx];
        N++;
      }
    }

    if(N) pVD[myID] /= (double) N; 

    memset(pVDTmp,0,sizeof(int)*iCells);
  }
  
  free(pTmp);
  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
  ENDVERBATIM
}









FUNCTION GetLoopLength () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetLoopLength ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetLoopLength ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetLoopLength ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);

  double* pLoop,  
          *pThru;   

  if( vector_arg_px(3,&pLoop) < iCells || vector_arg_px(4,&pThru) < iCells){
    printf("GetLoopLength ERRF: arg 3,4 must be Vectors with size >= %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }
  double dSubsamp = ifarg(5)?*getarg(5):1.0;

  unsigned int iSeed = ifarg(6)?(unsigned int)*getarg(6):INT_MAX-109754;

  
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetLoopLength ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }

  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0, iTmpSz = 0, jdx = 0;

  double  dgzt = 0.0; 
  int* pTmp = 0 , found = 0;
  double* pUse = 0; 
  
  if(dSubsamp<1.0){ 
     pUse = (double*)malloc(iCells*sizeof(double));
     mcell_ran4(&iSeed, pUse, iCells, 1.0);
  }

  pTmp = (int*)calloc(iCells,sizeof(int)); 

  if( verbose > 0 ) printf("searching loops from id: ");

  int* pVDTmp = (int*)calloc(iCells,sizeof(int)) , myID;

  for(myID=0;myID<iCells;myID++){

    if(!pLoop[myID]) continue;

    if(verbose > 0 && myID%1000==0)printf("%d ",myID); 

    
    if(pUse && pUse[myID]>=dSubsamp) continue;

    iCheckSz = 0; idx = 0; iDist = 1; youID = 0; youKidID = 0; found = 0;

    pVDTmp[myID]=1;

    
    for(idx=0;idx<pLen[myID];idx++){
      youID = pLV[myID][idx];
      if(youID==myID) {
        found = 1;
        pVD[myID]=iDist;
        iCheckSz=0;
        break;
      }
      if(pThru[youID] && !pVDTmp[youID]){
        pVDTmp[youID]=iDist;
        pCheck[iCheckSz++]=youID;
      }
    }

    iTmpSz = 0;  jdx=0;

    iDist++;
  
    
    while(iCheckSz>0){
      iTmpSz = 0;
      for(idx=0;idx<iCheckSz;idx++){
        youID=pCheck[idx];
        for(jdx=0;jdx<pLen[youID];jdx++){
          youKidID=pLV[youID][jdx];
          if(youKidID==myID){
            pVD[myID]=iDist;
            found = 1;
            break;
          }
          if(pThru[youKidID] && !pVDTmp[youKidID]){ 
            pTmp[iTmpSz++] = youKidID; 
            pVDTmp[youKidID]=iDist;
          }
        }
      }
      if(found) break;
      iCheckSz = iTmpSz;
      if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
      iDist++;
    }
    memset(pVDTmp,0,sizeof(int)*iCells);
  }
  
  free(pTmp);
  if(pUse) free(pUse); 
  free(pCheck);
  FreeListVec(&pList);  
  free(pVDTmp);

  if( verbose > 0 ) printf("\n");

  return 1.0;
  ENDVERBATIM
}







FUNCTION GetPathEV () {
  VERBATIM
  ListVec* pList = AllocListVec(*hoc_objgetarg(1));
  if(!pList){
    printf("GetPathEV ERRA: problem initializing first arg!\n");
    return 0.0;
  }
 
  int iCells = pList->isz; 
  if(iCells < 2){
    printf("GetPathEV ERRB: size of List < 2 !\n");
    FreeListVec(&pList);
    return 0.0;
  }

  double** pLV = pList->pv;
  int* pLen = pList->plen;

  
  double* pVD; 
  int iVecSz = vector_arg_px(2,&pVD) , i = 0;
  if(!pVD || iVecSz < iCells){
    printf("GetPathEV ERRE: arg 2 must be a Vector with size %d\n",iCells);
    FreeListVec(&pList);
    return 0.0;
  }  
  memset(pVD,0,sizeof(double)*iVecSz);

  
  int myID = (int) *getarg(3);
  if(myID < 0 || myID >= iCells){
    printf("GetPathEV ERRF: invalid id = %d\n",myID);
    FreeListVec(&pList);
    return 0.0;
  }

  
  int iStartID = ifarg(4) ? (int)*getarg(4) : 0,
      iEndID = ifarg(5) ? (int)*getarg(5) : iCells - 1,
      iMaxDist = ifarg(6)? (int)*getarg(6): -1;

  if(iStartID < 0 || iStartID >= iCells ||
     iEndID < 0 || iEndID >= iCells ||
     iStartID >= iEndID){
       printf("GetPathEV ERRH: invalid ids start=%d end=%d numcells=%d\n",iStartID,iEndID,iCells);
       FreeListVec(&pList);
       return 0.0;
     }

  
  if(iMaxDist==0){
    printf("GetPathEV ERRI: invalid maxdist=%d\n",iMaxDist);
    FreeListVec(&pList);
    return 0.0;
  }

  
  int* pCheck = (int*)malloc(sizeof(int)*iCells);
  if(!pCheck){
    printf("GetPathEV ERRG: out of memory!\n");
    FreeListVec(&pList);
    return 0.0;
  }
  int iCheckSz = 0, idx = 0, iDist = 1 , youID = 0, youKidID=0;

  pVD[myID]=1;

  
  for(idx=0;idx<pLen[myID];idx++){
    youID = pLV[myID][idx];
    if(youID>=iStartID && youID<=iEndID && !pVD[youID]){
      pVD[youID]=(double)iDist;
      pCheck[iCheckSz++]=youID;
    }
  }

  int* pTmp = (int*)malloc(sizeof(int)*iCells);
  int iTmpSz = 0 , jdx=0;

  iDist++;
  
  
  while(iCheckSz>0 && (iMaxDist==-1 || iDist<=iMaxDist)){
    iTmpSz = 0;
    for(idx=0;idx<iCheckSz;idx++){
      youID=pCheck[idx];
      for(jdx=0;jdx<pLen[youID];jdx++){
        youKidID=pLV[youID][jdx];
        if(youKidID >= iStartID && youKidID <=iEndID && !pVD[youKidID]){ 
          pTmp[iTmpSz++] = youKidID; 
          pVD[youKidID]=(double)iDist;
        }
      }
    }
    iCheckSz = iTmpSz;
    if(iCheckSz) memcpy(pCheck,pTmp,sizeof(int)*iCheckSz);
    iDist++;
  }

  pVD[myID]=0.0;
 
  free(pCheck);
  free(pTmp);
  FreeListVec(&pList);

  return 1.0;
  ENDVERBATIM
}


FUNCTION Factorial () {
  VERBATIM
  double N = (int)*getarg(1) , i = 0.0;
  double val = 1.0;
  if(N<=1) return 1.0;
  if(N>=171){
    double PI=3.1415926535897932384626433832795;
    double E=2.71828183;
    val=sqrt(2*PI*N)*(pow(N,N)/pow(E,N));
  } else {
    for(i=2.0;i<=N;i++) val*=i;
  }
  return (double) val;  
  ENDVERBATIM
}



FUNCTION perm () {
  VERBATIM
  if(ifarg(3)){
    double N = (int)*getarg(1);
    double R = (int)*getarg(2);
    double b = *getarg(3);
    double val = N/b;
    int i = 0;
    for(i=1;i<R;i++){
      N--;
      val*=(N/b);
    }
    return val;
  } else {
    int N = (int)*getarg(1);
    int R = (int)*getarg(2);
    int val = N;
    int i = 0;
    for(i=1;i<R;i++){
      N--;
      val*=N;
    }
    return (double)val;
  }
  ENDVERBATIM
}


PROCEDURE install () {
 if(INSTALLED==1){
   printf("Already installed $Id
 } else {
 INSTALLED=1
 VERBATIM
 install_vector_method("gzmean" ,gzmean);
 install_vector_method("nnmean" ,nnmean);
 install_vector_method("copynz" ,copynz);
 ENDVERBATIM
 printf("Installed $Id
 }
}