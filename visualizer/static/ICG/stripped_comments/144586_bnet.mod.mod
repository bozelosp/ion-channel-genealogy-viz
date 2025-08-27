NEURON {
  ARTIFICIAL_CELL BNET
  
  RANGE tstep
  RANGE xloc, yloc, zloc
  GLOBAL verbose, installed
  POINTER sop 
}

PARAMETER {
  tstep = 0
  verbose = 0
  sop = 0
}

ASSIGNED {
  installed
}



CONSTRUCTOR {
  VERBATIM 
  boonet* pnet;
  int sz; 
  if((sz = (int)*getarg(1)) < 1) {
    printf("BNET err0: must have a network with positive # of nodes!\n");
    hxe();
  }
  _p_sop = (void*) makeboonet(sz);
  pnet = BP;
  pnet->id = ifarg(2) ? (int) *getarg(2) : 0;
  ENDVERBATIM
}




FUNCTION setrule () {
  VERBATIM
  double *psrc, weight, *psrcstate;
  int targid, i, nsrc;
  if (! ifarg(4) ) {
    printf("BNET.setrule(vsource, targid, weight, vsourcestates) - sets a rule.\n");
    printf("vsource has node IDs of sources, targid is ID of target, weight is (-1,1) for inhib/activating\n");
    printf("vsourcestate has states the source must be in for the given rule to be turned on.\n");
    return 0.0;
  }
  if ( (nsrc = vector_arg_px(1,&psrc)) < 1) {
    printf("BNET.setrule WARN0: empty source Vector!\n");
    return 0.0;
  }
  targid = (int) *getarg(2);
  if(targid < 0 || targid >= BP->numnodes) {
    printf("BNET.setrule ERR0: invalid target id : %d\n",targid);
    return 0.0;
  }
  weight = (int) *getarg(3);
  if( nsrc != vector_arg_px(4,&psrcstate) ) {
    printf("BNET.setrule ERR1: vsource, vsourcestate must have same size!\n");
    return 0.0;
  }
  for(i=0;i<nsrc;i++) {
    if( psrc[i] < 0 || psrc[i] >= BP->numnodes) {
      printf("BNET.setrule ERR2: invalid source node id %d. netsize=%d\n",(int)psrc[i],BP->numnodes);
      return 0.0;
    }
    if(verbose>1) printf("adding rule from %d -> %d : w = %g\n",(int)psrc[i],targid,weight);
  }
  addrule(BP, psrc, targid, weight, psrcstate, nsrc);
  return 1.0;
  ENDVERBATIM
}


PROCEDURE clearrules () {
  VERBATIM
  BP->nrules = 0;
  ENDVERBATIM
}


PROCEDURE pr () {
  VERBATIM
  int i, j, k;
  bnode* pnodes = BP->pnodes;
  brule* prule;
  char srcstr[4096], stmp[4096];
  printf("net: numnodes=%d, numrules=%d\n",BP->numnodes,BP->nrules);
  for(i=0;i<BP->numnodes;i++) {
    if(pnodes[i].name) {
      printf("%s: state=%d, count=%d\n",pnodes[i].name,pnodes[i].state,pnodes[i].count);
    } else {
      printf("%d: state=%d, count=%d\n",i,pnodes[i].state,pnodes[i].count);
    }    
  }
  for(i=0;i<BP->nrules;i++) {
    prule = &BP->prules[i];
    srcstr[0]=0;
    for(j=0;j<prule->nsrc;j++) {
      if(prule->psrc[j]->name) {
        sprintf(stmp,"%s%s%s ", j>0?"AND ":"", prule->psrcstate[j]?"":"!", prule->psrc[j]->name);
      } else {
        sprintf(stmp,"%s%s%s ", j>0?"AND ":"", prule->psrcstate[j]?"":"!", prule->psrc[j]->id);
      }
      strcat(srcstr,stmp);
    }
    if(prule->ptarg->name) {
      printf("%s-> %s , w = %d\n",srcstr,prule->ptarg->name,prule->weight);
    } else {
      printf("%s-> %d , w = %d]\n",srcstr,prule->ptarg->id,prule->weight);
    }
  }
  ENDVERBATIM
}




FUNCTION graphviz () {
  VERBATIM
  int i, j, k, LR, fsz, w, h;
  bnode* pnodes = BP->pnodes;
  brule* prule;
  char *ncolor, *fcolor, *arrowtype, *lstyle, *shape;
  char buf[4096], *dotname, *fname, *ext, fontsize[128];
  double penw; 
  FILE* fp = 0x0;
  dotname = ifarg(1) ? gargstr(1) : 0x0;
  fname = ifarg(2) ? gargstr(2) : 0x0;
  ext =   ifarg(3) ? gargstr(3) : "gif";
  LR = ifarg(4) ? (int) *getarg(4) : 1;
  w = ifarg(5) ? (int) *getarg(5) : -1;
  h = ifarg(6) ? (int) *getarg(6) : -1;
  fsz = ifarg(7) ? (int) *getarg(7) : -1;
  if(fsz==-1) sprintf(fontsize,"%s"," "); else sprintf(fontsize,"fontsize=%d,",fsz);
  if(fname) if( !(fp = fopen(dotname,"w"))) {
    printf("BNET.graphviz ERR0: could not open %s\n",fname);
    return 0.0;
  }
  sprintf(buf, "%s", "digraph G {\n"); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
  if(LR){sprintf(buf, "%s", "\trankdir=LR;\n"); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); }
  if(w>0 && h>0) {sprintf(buf, "size=\"%d,%d\"\n",w,h); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf);}
  for(i=0;i<BP->numnodes;i++) {
    ncolor = BP->pnodes[i].knockout ? "white" : BP->pnodes[i].state > 0 ? "black" : "gray";
    fcolor = BP->pnodes[i].knockout ? "black" : "white";
    shape = pnodes[i].sthresh > 0 ? "invtriangle" : "doublecircle";
    if(BP->pnodes[i].name) {
      sprintf(buf,"\t%s [fontcolor=%s,%sstyle=filled,shape=%s,fillcolor=%s,color=%s]\n",
              BP->pnodes[i].name,fcolor,fontsize,shape,ncolor,ncolor);
    } else {
      sprintf(buf,"\t%d [fontcolor=%s,%sstyle=filled,shape=%s,fillcolor=%s,color=%s]\n",
              i,fcolor,fontsize,shape,ncolor,ncolor);
    }
    if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
  }
  for(i=0;i<BP->nrules;i++) {
    prule = &BP->prules[i];
    for(j=0;j<prule->nsrc;j++) {
      penw = prule->psrcstate[j] == prule->psrc[j]->state ? 6.0 : 1.0;
      arrowtype = prule->weight < 0 ? "tee" : "open";
      lstyle = prule->psrcstate[j] == 0 ? ",style=dashed" : " ";
      if(prule->psrc[j]->name) {
        if(prule->ptarg->name) {
          sprintf(buf,"\t%s -> %s [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->name,prule->ptarg->name,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
        } else {
          sprintf(buf,"\t%s -> %d [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->name,prule->ptarg->id,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
        }
      } else if(prule->ptarg->name) {
          sprintf(buf,"\t%d -> %s [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->id,prule->ptarg->name,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
      } else {
          sprintf(buf,"\t%d -> %d [arrowhead=%s,penwidth=%g,color=%s%s]\n",
                  prule->psrc[j]->id,prule->ptarg->id,arrowtype,penw,prule->weight>0?"red":"blue",lstyle);
      }
      if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
    }
  }
  sprintf(buf,"%s","}\n"); if(fp) fprintf(fp,"%s",buf); else fprintf(stdout,"%s",buf); 
  if(fp) fclose(fp);
  if(fname) {
    sprintf(buf,"dot %s -T%s > %s",dotname,ext,fname);
    if(0!=system(buf)) {printf("BNET.graphviz ERR1 : couldn't run %s\n",buf); return 0.0;}
  }
  return 1.0;
  ENDVERBATIM
}

DESTRUCTOR {
  VERBATIM
  freeboonet(BP);
  ENDVERBATIM
}


PROCEDURE start () {
  tstep = 0
  VERBATIM
  startboonet(BP); 
  ENDVERBATIM
}


INITIAL {
  start()
}


PROCEDURE strvalid () {
  VERBATIM
  char *pname;
  static char *pnames[6] = {"state", "count", "knockout", "start", "sthresh", "scount"};
  int i;
  pname = gargstr(1);
  for(i=0;i<6;i++) {
    if(!strcmp(pname,pnames[i])) return i;
  }
  return -1;
  ENDVERBATIM
}



PROCEDURE getscount () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getscount(vec) - returns scount of each node in vec\n");
    return;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].scount;
  ENDVERBATIM
}



FUNCTION setscount () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setscount(vec) - sets scount of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setscount ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].scount = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}



PROCEDURE getsthresh () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getsthresh(vec) - returns sthresh of each node in vec\n");
    return;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].sthresh;
  ENDVERBATIM
}



FUNCTION setsthresh () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setsthresh(vec) - sets sthresh of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setsthresh ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].sthresh = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}



PROCEDURE getcount () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getcount(vec) - returns count of each node in vec\n");
    return;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].count;
  ENDVERBATIM
}



FUNCTION setcount () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setcount(vec) - sets count of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setcount ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].count = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}



PROCEDURE getstate () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getstate(vec) - returns state of each node in vec\n");
    return;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].state;
  ENDVERBATIM
}



FUNCTION setstate () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setstate(vec) - sets state of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setstate ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].state = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}



PROCEDURE getstart () {
  VERBATIM
  double *ps; int i; void *vs;
  if(!ifarg(1)) {
    printf("BNET.getstart(vec) - returns start state of each node in vec\n");
    return;
  }
  vs = vector_arg(1);
  ps = vector_newsize(vs,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) ps[i] = (double) BP->pnodes[i].start;
  ENDVERBATIM
}



FUNCTION setstart () {
  VERBATIM
  double *ps; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.setstart(vec) - sets start state of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&ps)) != BP->numnodes ) {
    printf("BNET.setstart ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].start = (int) ps[i];
  return 1.0;
  ENDVERBATIM
}



FUNCTION setknockout () {
  VERBATIM
  double *pk; int i, sz;
  if(!ifarg(1)) {
    printf("BNET.knockout(vec) - sets knockout flag of each node in vec\n");
    return 0.0;
  }
  if( (sz = vector_arg_px(1,&pk)) != BP->numnodes ) {
    printf("BNET.knockout ERR0: vec.size(%d) != BNET.numnodes(%d)\n",sz,BP->numnodes);
    return 0.0;
  }
  for(i=0;i<BP->numnodes;i++) BP->pnodes[i].knockout = (int) pk[i];
  return 1.0;
  ENDVERBATIM
}



PROCEDURE getknockout () {
  VERBATIM
  double *pk; int i; void *vk;
  if(!ifarg(1)) {
    printf("BNET.getknockout(vec) - returns knockout flag of each node in vec\n");
    return;
  }
  vk = vector_arg(1);
  pk = vector_newsize(vk,BP->numnodes);
  for(i=0;i<BP->numnodes;i++) pk[i] = (double) BP->pnodes[i].knockout;
  ENDVERBATIM
}


FUNCTION setnname () {
  VERBATIM
  int id, sz; char *name;
  id = (int) *getarg(1);
  if(id < 0 || id >= BP->numnodes) {
    printf("BNET.setnname ERR0: invalid node index %d\n",id);
    return 0.0;    
  }
  name = gargstr(2);
  if(!(sz=strlen(name))) {
    printf("BNET.setnname ERR1: empty string\n");
    return 0.0;
  }
  if(BP->pnodes[id].name) free(BP->pnodes[id].name);
  if(!(BP->pnodes[id].name = (char*) malloc(sizeof(char) * (sz + 1)))) {
    printf("BNET.setnname ERR2: couldn't alloc mem for %s\n",name);
    return 0.0;
  }
  strcpy(BP->pnodes[id].name,name);
  return 1.0;
  ENDVERBATIM
}


FUNCTION getnname () {
  VERBATIM
  int i, id, sz; char **pname, string[BUFSIZ];
  char** hoc_pgargstr();
  id = (int) *getarg(1);
  if(id < 0 || id >= BP->numnodes) {
    printf("BNET.getnname ERR0: invalid node index %d\n",id);
    return 0.0;    
  }
  if(!BP->pnodes[id].name || !(sz=strlen(BP->pnodes[id].name))) {
    printf("BNET.getnname ERR1: node %d has no name\n",id);
    return 0.0;
  }
  for(i=0;i<sz && i<BUFSIZ;i++) string[i] = BP->pnodes[id].name[i];
  if(i < BUFSIZ) string[i]=0;
  printf("Aname is %s, %s\n",BP->pnodes[id].name,string);
  pname = hoc_pgargstr(2);
  printf("Bname is %s, %s\n",BP->pnodes[id].name,string);
  hoc_assign_str(pname,string);
  return 1.0;
  ENDVERBATIM
}


FUNCTION advancebn () {
  VERBATIM
  advanceboonet(BP);  
  tstep = tstep + 1;
  return tstep;
  ENDVERBATIM
}

FUNCTION advancebnfor () {
  VERBATIM
  int i, n;
  n = (int) *getarg(1);
  for(i=0;i<n;i++) advancebn();
  return tstep;
  ENDVERBATIM
}

FUNCTION numnodes () {
  VERBATIM
  return BP->numnodes;
  ENDVERBATIM
}

FUNCTION numrules () {
  VERBATIM
  return BP->nrules;
  ENDVERBATIM
}


FUNCTION id () {
  VERBATIM
  return (double) BP->id;
  ENDVERBATIM
}