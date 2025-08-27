NEURON {
  ARTIFICIAL_CELL JitCon
  POINTER sop                          
}

PARAMETER {
  sop=0
}

ASSIGNED {
}



CONSTRUCTOR {
  VERBATIM 
  { int lid,lty,lco,i;
    if (ifarg(1)) { lid=(int) *getarg(1); } else { lid= UINT_MAX; }
    if (ifarg(2)) { lty=(int) *getarg(2); } else { lty= -1; }
    if (ifarg(3)) { lco=(int) *getarg(3); } else { lco= -1; }
    _p_sop = (void*)ecalloc(1, sizeof(id0)); 
    ip = IDP;
    ip->id=lid; ip->type=lty; ip->col=lco; ip->pg=0x0; ip->dvi=0x0; ip->sprob=0x0;
    ip->jcn = 0;
    process=(int)getpid();
    CNAME[SU]="SU"; CNAME[DP]="DP"; CNAME[IN]="IN";
  }
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM { 
  free(IDP);
  }
  ENDVERBATIM
}


INITIAL {
}


NET_RECEIVE (w) { LOCAL tmp,jcn
  VERBATIM
  ip = IDP;
  _ljcn=ip->jcn;
  tpnt = _pnt; 
  
  
  if (! ip->jcn) { 
    net_event(tpnt, t); 
  
  } else if (_lflag==0) { 
    jitcon(t);
  
  } else if (_lflag<0) { 
    callback(_lflag); 
  }
  ENDVERBATIM
}








PROCEDURE jitcon (tm) {
  VERBATIM {
  double mindel, randel, idty, *x; int prty, poty, i, j, k, dv, dvt; 
  Point_process *pnt; void* voi;
  
  
  ip=IDP;
  prty=(int)ip->type;
  if (ip->jcn==1) { 
    if (!pg) {printf("No network defined -- must run jitcondiv()\n"); hxe();}
#if defined(t)
    if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,ip->del[0]+t,-1.); 
#else
    if (ip->dvt>0) net_send((void**)0x0, wts,tpnt,ip->del[0],-1.); 
#endif
  } else if (ip->jcn==3) { mkdvi(); 
  } else if (ip->jcn==2) { 
    
    
    
    
    
    for (i=0,k=0,dvt=0;i<CTYN;i++) { 
      poty=cty[i];
      printf("%d:%d:%d ",prty,poty,(int)DVG(prty,poty));
      dvt+=(int)DVG(prty,poty);
    }
    for (i=0,k=0,dvt=0;i<CTYN;i++) { 
      if (dv>0) {
        mcell_ran4(&sead, &randel , 1, 2.*DELD(prty,poty)); 
        randel+=((_ltm-t)+DELM(prty,poty)-DELD(prty,poty));
        if (randel<=0) randel= -randel;
        if (dv>scrsz) {
          printf("A:Divergence exceeds scrsz: %d>%d for %d->%d\n",dv,scrsz,prty,poty); hxe(); }
        mcell_ran4(&sead, scr ,  dv, pg->ixe[poty]-pg->ix[poty]+1); 
        if (ifarg(2)) { 
          voi=vector_arg(2);
          x=vector_newsize(voi,k+dv+1);
          x[k++]= -randel; 
          for (j=0;j<dv;j++,k++) x[k]=floor(scr[j]+pg->ix[poty]);
        } else for (j=0;j<dv;j++) {
          pnt=(Point_process *)(ivoc_list_item(ce, (int)(scr[j]+pg->ix[poty])))->u.this_pointer;
          idty=(double)(FOFFSET+ip->id)+0.1*(double)ip->type+0.01;
#if defined(t)
          net_send((void**)0x0, wts, pnt, randel+t, idty);
#else
          net_send((void**)0x0, wts, pnt, randel, idty);
#endif
        }
      }
    }
  } 
  }
  ENDVERBATIM  
}

PROCEDURE callback (fl) {
  VERBATIM {
  int ii,jj; double idty, del; double weed, prid, prty, poid, poty, w;
  unsigned int valseed;
  ii=(unsigned int)((-_lfl)-1); 
  ip=IDP;
  idty=(double)(FOFFSET+ip->id)+0.1*(double)ip->type+0.01;
  jj=ii+1;
  if (jj<ip->dvt) {
    del= ip->del[jj] - ip->del[ii];
#if defined(t)
    net_send((void**)0x0, wts,tpnt,del+t,(double) -(jj+1)); 
#else
    net_send((void**)0x0, wts,tpnt,del,(double) -(jj+1)); 
#endif
  }
  if (ip->sprob[ii]) { 
    
    
    poid=ip->dvi[ii]->_prop->param[3]; 
    poty=ip->dvi[ii]->_prop->param[4];
    prid=ip->dvi[ii]->_prop->param[5];
    prty=ip->dvi[ii]->_prop->param[6];
    if ((((double)ip->type)!=prty) || (((double)ip->id)!=prid)) {
      printf("callback() ERR: %g %g %g %g %g %g\n",\
             ((double)ip->type),prty,((double)ip->id),prid,poty,poid); hxe(); }
    weed=prty*allcells+poty*100+prid*10+poid;
    
    valseed=(unsigned int)weed; w=WMAT((int)prty,(int)poty);
    mcell_ran4(&valseed, &wts, 1, 2*0.01*w); 
    
    wts[0]+=0.99*w; 
                    
    (*pnt_receive[ip->dvi[ii]->_prop->_type])(ip->dvi[ii], wts, 0.0); } 
  } 
  ENDVERBATIM
}


PROCEDURE mkdvi () {
VERBATIM {
  int i,j,k,prty,poty,dv,dvt,dvii; double *x, *db, *dbs; 
  Object *lb;  Point_process *pnnt, **da, **das;
  ip=IDP; ip->pg=pg; 
  prty=ip->type;
  sead=((unsigned int)ip->id)*1e6;
  for (i=0,k=0,dvt=0;i<CTYN;i++) { 
    poty=cty[i];
    dvt+=DVG(prty,poty);
  }
  da =(Point_process **)malloc(dvt*sizeof(Point_process *));
  das=(Point_process **)malloc(dvt*sizeof(Point_process *)); 
  db =(double *)malloc(dvt*sizeof(double)); 
  dbs=(double *)malloc(dvt*sizeof(double)); 
  for (i=0,k=0,dvii=0;i<CTYN;i++) { 
    poty=cty[i];
    dv=DVG(prty,poty);
    if (dv>0) {
      sead+=dv;
      if (dv>scrsz) {
        printf("B:Divergence exceeds scrsz: %d>%d for %d->%d\n",dv,scrsz,prty,poty); hxe(); }
      mcell_ran4(&sead, scr ,  dv, pg->ixe[poty]-pg->ix[poty]+1);
      for (j=0;j<dv;j++) {
        if (!(lb=ivoc_list_item(ce,(unsigned int)floor(scr[j]+pg->ix[poty])))) {
          printf("JitCon:callback %g exceeds %d for list ce\n",floor(scr[j]+pg->ix[poty]),cesz); 
          hxe(); }
        pnnt=(Point_process *)lb->u.this_pointer;
        da[j+dvii]=pnnt;
      }
      mcell_ran4(&sead, scr , dv, 2*DELD(prty,poty));
      for (j=0;j<dv;j++) {
        db[j+dvii]=scr[j]+DELM(prty,poty)-DELD(prty,poty); 
        if (db[j+dvii]<0) db[j+dvii]=-db[j+dvii];
      }
      dvii+=dv;
    }
  }
  gsort2(db,da,dvt,dbs,das);
  ip->del=dbs;   ip->dvi=das;   ip->dvt=dvt;
  ip->sprob=(unsigned char *)malloc(dvt*sizeof(char *)); 
  for (i=0;i<dvt;i++) ip->sprob[i]=1; 
  free(da); free(db); 
  }
ENDVERBATIM
}


FUNCTION getdvi () {
  VERBATIM 
  {
  int j,dvt; double *dbs, *x;
  void* voi; Point_process **das;
  ip=IDP; ip->pg=pg; 
  dvt=ip->dvt;
  dbs=ip->del;   das=ip->dvi;
  _lgetdvi=(double)dvt; 
  if (!ifarg(1)) return _lgetdvi; 
  if (hoc_is_double_arg(1)) { 
    
  }
  voi=vector_arg(1);
  x=vector_newsize(voi,dvt);
  for (j=0;j<dvt;j++) {
    x[j]=(double)das[j]->_prop->param[3]; 
  }
  voi=vector_arg(2);
  x=vector_newsize(voi,dvt);
  for (j=0;j<dvt;j++) x[j]=dbs[j];
  if (ifarg(3)) {
    voi=vector_arg(3);
    x=vector_newsize(voi,dvt);
    for (j=0;j<dvt;j++) x[j]=(double)ip->sprob[j];
  }
  }
ENDVERBATIM
}


PROCEDURE setdvi () {
VERBATIM {
  int i,j,k,dvt,dvu,ddvi,lbcnt; double *y, *db, *dbs, id;
  void* voi; Object *lb, *lc, *syo; Point_process *pnnt, **da, **das;
  ip=IDP; ip->pg=pg; id=(double)ip->id;
  dvt=vector_arg_px(1, &y); 
  i=vector_arg_px(2, &db);
  if (i != dvt) {printf("setdvi() ERR vec sizes: %d %d\n",dvt,j); hxe();}
  if (ip->type==IN) syo=gal; else syo=aml; 
  da =(Point_process **)malloc(dvt*sizeof(Point_process *));
  das=(Point_process **)malloc(dvt*sizeof(Point_process *)); 
  dbs=(double *)malloc(dvt*sizeof(double)); 
  for (j=0,i=0;j<dvt;j++) {
    lb=ivoc_list_item(syo,(unsigned int)y[j]); 
    if (!lb) { printf("JitCon:callback %g exceeds %d for list ce\n",y[j],cesz); hxe(); }
    lbcnt=ivoc_list_count(lb);
    for (k=0;k<lbcnt;k++) {
      lc=ivoc_list_item(lb,k);
      if (((Point_process *)lc->u.this_pointer)->_prop->param[5] == id) {
        da[j]=(Point_process *)lc->u.this_pointer;
        break;
      }
    }
    if (k==lbcnt) { printf("setdvi() ERR: id %g not found\n",id); hxe(); }
  }
  gsort2(db,da,dvt,dbs,das);
  ip->del=dbs;   ip->dvi=das;   ip->dvt=dvt;
  ip->sprob=(unsigned char *)malloc(dvt*sizeof(char *)); 
  for (j=0;j<dvt;j++) ip->sprob[j]=1; 
  free(da);
  }
ENDVERBATIM
}



PROCEDURE prune () {
  VERBATIM 
  {
  double *x, p; int nx,j;
  ip=IDP; ip->pg=pg;
  if (hoc_is_double_arg(1)) {
    p=*getarg(1);
    if (p<0 || p>1) {printf("JitCon:pruneERR0:need # [0,1] to prune [ALL,NONE]: %g\n",p); hxe();}
    if (p==1.) printf("JitConpruneWARNING: pruning 100% of cell %d\n",ip->id);
    if (ip->dvt>scrsz) {
      printf("JitConpruneB:Div exceeds scrsz: %d>%d\n",ip->dvt,scrsz); hxe(); }
    for (j=0;j<ip->dvt;j++) ip->sprob[j]=1; 
    if (p==0.) return; 
    sead=(ifarg(2))?(unsigned int)*getarg(2):(unsigned int)ip->id*1e6;
    mcell_ran4(&sead, scr , ip->dvt, 1.0); 
    for (j=0;j<ip->dvt;j++) if (scr[j]<p) ip->sprob[j]=0; 
  } else {
    nx=vector_arg_px(1,&x);
    if (nx!=ip->dvt) {printf("JitCon:pruneERRA:Wrong size vector:%d!=%d\n",nx,ip->dvt); hxe();}
    for (j=0;j<ip->dvt;j++) ip->sprob[j]=(unsigned char)x[j];
  }
  }
ENDVERBATIM
}

VERBATIM 

int gsort2 (double *db, Point_process **da,int dvt,double *dbs, Point_process **das) {
  unsigned int *scr, i;
  scr=scrset(dvt);
  for (i=0;i<dvt;i++) scr[i]=i;
  nrn_mlh_gsort(db, scr, dvt, cmpdfn);
  for (i=0;i<dvt;i++) {
    dbs[i]=db[scr[i]];
    das[i]=da[scr[i]];
  }
}
ENDVERBATIM

PROCEDURE freedvi () {
  VERBATIM
  { 
    int i, poty; id0 *jp;
    jp=IDP;
    if (jp->dvi) {
      free(jp->dvi);
      free(jp->del);
      jp->dvi=0x0;
      jp->del=0x0;
    }
  }
  ENDVERBATIM
}

FUNCTION qstats () {
  VERBATIM {
    double stt[3]; int lct,flag; FILE* tf;
    if (ifarg(1)) {tf=hoc_obj_file_arg(1); flag=1;} else flag=0;
    lct=cty[IDP->type];
    _lqstats = nrn_event_queue_stats(&stt);
    printf("QUEUE: Inserted %g; removed %g\n",stt[0],stt[2]);
    if (flag) {
      fprintf(tf,"QUEUE: Inserted %g; removed %g remaining: %g\n",stt[0],stt[2],_lqstats);
    }
  }
  ENDVERBATIM
}

FUNCTION qsz () {
  VERBATIM {
    double stt[3];
    _lqsz = nrn_event_queue_stats(&stt);
  }
  ENDVERBATIM
}

PROCEDURE qclr () {
  VERBATIM {
    clear_event_queue();
  }
  ENDVERBATIM
}


PROCEDURE jitcondiv () {
  VERBATIM {
  Symbol *sym; int i,j; char name[100];
  ip->pg=(postgrp *)malloc(sizeof(postgrp));
  pg=ip->pg;
  sym = hoc_lookup("ce"); ce = (*(hoc_objectdata[sym->u.oboff].pobj));
  sym = hoc_lookup("aml"); aml = (*(hoc_objectdata[sym->u.oboff].pobj));
  sym = hoc_lookup("gal"); gal = (*(hoc_objectdata[sym->u.oboff].pobj));
  cesz = ivoc_list_count(ce);
  if (cesz!=(i=ivoc_list_count(aml)) || cesz!=(j=ivoc_list_count(gal))) {
    printf("All 3 lists should be same size: ce,aml,gal %d,%d,%d\n",cesz,i,j); hxe(); }
  cty[0]=SU; cty[1]=IN; 
  CTYPi=HVAL("CTYPi"); STYPi=HVAL("STYPi"); scrsz=HVAL("scrsz"); allcells=HVAL("allcells");
  pg->ix =HPTR("ix"); pg->ixe=HPTR("ixe"); 
  pg->dvg=HPTR("div"); 
  pg->wmat=HPTR("wmat");
  pg->delm=HPTR("delm"); pg->deld=HPTR("deld");
  scr=HPTR("scr");
  if (!ce) {printf("JitCon jitcondiv ERRA: ce not found\n"); hxe();}
  
  printf("Checking for possible seg error in double arrays: CTYPi==%d: ",CTYPi);
  
  printf("%d %d %d ",DVG(CTYPi-1,CTYPi-1),(int)pg->ix[CTYPi-1],(int)pg->ixe[CTYPi-1]);
  printf("%g ",WMAT(CTYPi-1,CTYPi-1));
  printf("%g %g ",DELM(CTYPi-1,CTYPi-1),DELD(CTYPi-1,CTYPi-1));
  printf("%d %g\n",scrsz,scr[scrsz-1]); 
  }
  ENDVERBATIM  
}


PROCEDURE probejcd () {
  VERBATIM {  int i,a[4];
    for (i=1;i<=3;i++) a[i]=(int)*getarg(i);
    printf("CTYPi: %d, STYPi: %d, ",CTYPi,STYPi);
    
    printf("wmat: %g\n",WMAT(a[1],a[2]));
  }
  ENDVERBATIM  
}


PROCEDURE vers () {
  printf("$Id
}

VERBATIM



static double* lop (Object *ob, unsigned int i) {
  Object *lb;
  lb = ivoc_list_item(ob, i);
  if (! lb) { printf("JitCon:lop %d exceeds %d for list ce\n",i,cesz); hxe();}
  pmt=ob2pntproc(lb);
  qp=*((id0**) &((pmt->_prop->dparam)[2])); 
  return pmt->_prop->param;
}


int stoppo () {
}
ENDVERBATIM



PROCEDURE lof () {
VERBATIM {
  Object *ob; int num,i,ii,j,k,si,nx;  double *vvo[7], *par; void *vv[7];
  ob = *(hoc_objgetarg(1));
  si=(int)*getarg(2);
  num = ivoc_list_count(ob);
  if (num!=7) { printf("JitCon lof ERR %d>7\n",num); hxe(); }
  for (i=0;i<num;i++) { 
    j = list_vector_px3(ob, i, &vvo[i], &vv[i]);
    if (i==0) nx=j;
    if (j!=nx) { printf("JitCon lof ERR %d %d\n",j,nx); hxe(); }
  }
 }
ENDVERBATIM
}


PROCEDURE ldv () {
VERBATIM {
  Object *ob; 
  int dv,num,i,j,ii,prty,poty,nx,offset;  
  double *vvo[100000], *x; void *vv[100000];
  ob = *(hoc_objgetarg(1));
  nx = vector_arg_px(2, &x); 
  prty=(int)*getarg(3);
  poty=(int)*getarg(4);
  if (ifarg(5)) offset=(int)*getarg(5); else offset=0;
  num = ivoc_list_count(ob);
  if (num!=nx) {printf("JitCon ldv ERRD %d %d\n",num,nx); hxe(); }
  
  
  if (num>1e5) { printf("JitCon ldv ERRA %d>1e5\n",num); hxe(); }
  i=0; nx=list_vector_px3(ob, i, &vvo[i], &vv[i]);
  dv=DVG(prty,poty);
  if (nx!=dv) { printf("JitCon ldv ERRB %d %d\n",dv,nx); hxe(); }
  for (i=1;i<num;i++) { 
    j = list_vector_px3(ob, i, &vvo[i], &vv[i]);    
    if (j!=nx) { printf("JitCon ldv ERRC %d %d\n",j,nx); hxe(); }
  }
  if (ii!=num) printf("INF ldv WARNING: only filled %d of %d columns\n",ii,num);
 }
ENDVERBATIM
}



PROCEDURE chk (f) {
  VERBATIM 
  {int i,lfg;
  lfg=(int)_lf;
  ip=IDP;
  if (lfg==1) {
  }
  if (lfg==2) { 
  }
  if (lfg==3) { 
  }
  if (lfg==4) { 
  }
  if (lfg==5) { 
  }}
  ENDVERBATIM
}


FUNCTION pid () {
  VERBATIM 
  printf("JitCon%d(%d/%d@%g) ",IDP->id,IDP->type,IDP->col,t);
  _lpid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION id () {
  VERBATIM
  if (ifarg(1)) IDP->id = (unsigned int) *getarg(1);
  _lid = (double)IDP->id;
  ENDVERBATIM
}

FUNCTION type () {
  VERBATIM
  if (ifarg(1)) IDP->type = (unsigned char) *getarg(1);
  _ltype = (double)IDP->type;
  ENDVERBATIM
}

FUNCTION col () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->col = (unsigned char) *getarg(1);
  _lcol = (double)ip->col;
  ENDVERBATIM
}

FUNCTION dbx () {
  VERBATIM 
  ip = IDP; 
  if (ifarg(1)) ip->dbx = (unsigned char) *getarg(1);
  _ldbx = (double)ip->dbx;
  ENDVERBATIM
}


FUNCTION fflag () { fflag=1 }




FUNCTION flag () {
  VERBATIM 
  {char *sf; int ii,i; unsigned char val;
  ip = IDP;
  sf = gargstr(1);
  for (ii=0;ii<iflnum && strncmp(sf, &iflags[ii*4], 3)!=0;ii++) ;
  if (ii==iflnum) {printf("JitCon ERR: %s not found as a flag (%s)\n",sf,iflags); hxe();}
  if (ifarg(2)) (&ip->type)[ii] = val = (unsigned char) *getarg(2);  
  _lflag=(double)(unsigned char)(&ip->type)[ii];
  if (ifarg(3)) for (i=0;i<cesz;i++) { 
    lop(ce,i); 
    (&qp->type)[ii]=val;
  }
  }
  ENDVERBATIM
}