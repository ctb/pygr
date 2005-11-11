
#include "intervaldb.h"

IntervalMap *read_intervals(int n,FILE *ifile)
{
  int i=0,start,end,target_id,target_start,target_end;
  IntervalMap *im=NULL;
  CALLOC(im,n,IntervalMap); /* ALLOCATE THE WHOLE ARRAY */
  while (fscanf(ifile," %d %d %d %d %d",&im[i].start,&im[i].end,
		&im[i].target_id,&im[i].target_start,
		&im[i].target_end)==5) {
    im[i].sublist= -1; /* DEFAULT: NO SUBLIST */
    i++;
  }
  if (i!=n) {
    fprintf(stderr,"Error: number of records read %d does not match allocation %d\n",i,n);
    exit(1);
  }
  return im;
 handle_malloc_failure:
  return NULL; /* NO CLEANUP ACTIONS REQUIRED */
}







int im_qsort_cmp(const void *void_a,const void *void_b)
{
  int a_start,a_end,b_start,b_end;
  IntervalMap *a=(IntervalMap *)void_a,*b=(IntervalMap *)void_b;
  SET_INTERVAL_POSITIVE(*a,a_start,a_end);
  SET_INTERVAL_POSITIVE(*b,b_start,b_end);
  if (a_start<b_start)
    return -1;
  else if (a_start>b_start)
    return 1;
  else if (a_end>b_end) /* SAME START: PUT LONGER INTERVAL 1ST */
    return -1;
  else if (a_end<b_end) /* CONTAINED INTERVAL SHOULD FOLLOW LARGER INTERVAL*/
    return 1;
  else
    return 0;
}


int sublist_qsort_cmp(const void *void_a,const void *void_b)
{
  IntervalMap *a=(IntervalMap *)void_a,*b=(IntervalMap *)void_b;
  if (a->sublist<b->sublist)
    return -1;
  else if (a->sublist>b->sublist)
    return 1;
  else if (a->start < b->start)
    return -1;
  else if (a->start > b->start)
    return 1;
  else
    return 0;
}



SublistHeader *build_nested_list(IntervalMap im[],int n,
				 int *p_n,int *p_nlists)
{
  int i=0,j,k,l,parent,nsub=0,nlists=0;
  IntervalMap *imsub=NULL;
  SublistHeader *subheader=NULL;

  qsort(im,n,sizeof(IntervalMap),im_qsort_cmp); /* SORT BY start, CONTAINMENT */
  while (i<n) { /* TOP LEVEL LIST SCAN */
    parent=i;
    i=parent+1;
    while (i<n && parent>=0) { /* RECURSIVE ALGORITHM OF ALEX ALEKSEYENKO */
      if (END_POSITIVE(im[i])<=END_POSITIVE(im[parent])) {/* i CONTAINED IN parent*/
	im[i].sublist=parent; /* MARK AS CONTAINED IN parent */
	nsub++; /* COUNT TOTAL #SUBLIST ENTRIES */
	parent=i; /* AND PUSH ONTO RECURSIVE STACK */
	i++; /* ADVANCE TO NEXT INTERVAL */
      }
      else /* POP RECURSIVE STACK*/
	parent=im[parent].sublist;
    }
  } /* AT THIS POINT sublist IS EITHER -1 IF NOT IN SUBLIST, OR INDICATES parent*/

  if (nsub>0) { /* WE HAVE SUBLISTS TO PROCESS */
    CALLOC(imsub,nsub,IntervalMap); /* TEMPORARY ARRAY FOR REPACKING SUBLISTS */
    for (i=j=0;i<n;i++) { /* GENERATE LIST FOR SORTING; ASSIGN HEADER INDEXES*/
      parent=im[i].sublist;
      if (parent>=0)  {/* IN A SUBLIST */
	imsub[j].start=i;
	imsub[j].sublist=parent;
	j++;
	if (im[parent].sublist<0) /* A NEW PARENT! SET HIS SUBLIST HEADER INDEX */
	  im[parent].sublist=nlists++;
      }
      im[i].sublist= -1; /* RESET TO DEFAULT VALUE: NO SUBLIST */
    }
    qsort(imsub,nsub,sizeof(IntervalMap),sublist_qsort_cmp);
    /* AT THIS POINT SUBLISTS ARE GROUPED TOGETHER, READY TO PACK */

    CALLOC(subheader,nlists,SublistHeader); /* SUBLIST HEADER INDEX */
    for (i=0;i<nsub;i++) { /* COPY SUBLIST ENTRIES TO imsub */
      j=imsub[i].start;
      parent=imsub[i].sublist;
      memcpy(imsub+i,im+j,sizeof(IntervalMap)); /* COPY INTERVAL */
      k=im[parent].sublist;
      if (subheader[k].len==0) /* START A NEW SUBLIST */
	subheader[k].start=i;
      subheader[k].len++; /* COUNT THE SUBLIST ENTRIES */
      im[j].start=im[j].end= -1; /* MARK FOR DELETION */
    } /* DONE COPYING ALL SUBLISTS TO imsub */

    for (i=j=0;i<n;i++) /* COMPRESS THE LIST TO REMOVE SUBLISTS */
      if (im[i].start!= -1 || im[i].end!= -1) { /* NOT IN A SUBLIST, SO KEEP */
	if (j<i) /* COPY TO NEW COMPACTED LOCATION */
	  memcpy(im+j,im+i,sizeof(IntervalMap));
	j++;
      }

    memcpy(im+j,imsub,nsub*sizeof(IntervalMap)); /* COPY THE SUBLISTS */
    for (i=0;i<nlists;i++) /* ADJUST start ADDRESSES FOR SHIFT*/
      subheader[i].start += j; 
    FREE(imsub);
	 *p_n = j; /* COPY THE COMPRESSED LIST SIZES BACK TO CALLER*/
  }
  else 
	 *p_n = n;
  *p_nlists=nlists;
  return subheader;
 handle_malloc_failure:
  FREE(imsub);  /* FREE ANY MALLOCS WE PERFORMED*/
  FREE(subheader);
  return NULL;
}


IntervalMap *interval_map_alloc(int n)
{
  IntervalMap *im=NULL;
  CALLOC(im,n,IntervalMap);
  return im;
 handle_malloc_failure:
  return NULL;
}



int find_overlap_start(int start,int end,IntervalMap im[],int n)
{
  int l=0,mid,r;

  r=n-1;
  while (l<r) {
    mid=(l+r)/2;
    if (END_POSITIVE(im[mid])<=start)
      l=mid+1;
    else
      r=mid;
  }
  if (l<n && HAS_OVERLAP_POSITIVE(im[l],start,end))
    return l; /* l IS START OF OVERLAP */
  else 
    return -1; /* NO OVERLAP FOUND */
}




int find_index_start(int start,int end,IntervalIndex im[],int n)
{
  int l=0,mid,r;

  r=n-1;
  while (l<r) {
    mid=(l+r)/2;
    if (END_POSITIVE(im[mid])<=start)
      l=mid+1;
    else
      r=mid;
  }
  return l; /* l IS START OF POSSIBLE OVERLAP */
}



int find_suboverlap_start(int start,int end,int isub,IntervalMap im[],
			  SublistHeader subheader[],int nlists)
{
  int i;

  if (isub>=0) {
    i=find_overlap_start(start,end,im+subheader[isub].start,subheader[isub].len);
    if (i>=0)
      return i+subheader[isub].start;
  }
  return -1;
}


IntervalIterator *interval_iterator_alloc()
{
  IntervalIterator *it=NULL;
  CALLOC(it,1,IntervalIterator);
  return it;
 handle_malloc_failure:
  return NULL;
}

int free_interval_iterator(IntervalIterator *it)
{
  IntervalIterator *it2,*it_next;
  FREE_ITERATOR_STACK(it,it2,it_next);
  return 0;
}


IntervalIterator *reset_interval_iterator(IntervalIterator *it)
{
  ITERATOR_STACK_TOP(it);
  it->n=0;
  return it;
}


void reorient_intervals(int n,IntervalMap im[],int ori_sign)
{
  int i,tmp;
  for (i=0;i<n;i++) {
    if ((im[i].start>=0 ? 1:-1)!=ori_sign) { /* ORIENTATION MISMATCH */
      tmp=im[i].start; /* SO REVERSE THIS INTERVAL MAPPING */
      im[i].start= -im[i].end;
      im[i].end =  -tmp;
      tmp=im[i].target_start;
      im[i].target_start= -im[i].target_end;
      im[i].target_end =  -tmp;
    }
  }
}

IntervalIterator *find_intervals(IntervalIterator *it0,int start,int end,
				 IntervalMap im[],int n,
				 SublistHeader subheader[],int nlists,
				 IntervalMap buf[],int nbuf,
				 int *p_nreturn)
{
  IntervalIterator *it=NULL,*it2=NULL;
  int ibuf=0,j,k,ori_sign=1;
  if (!it0) { /* ALLOCATE AN ITERATOR IF NOT SUPPLIED*/
    CALLOC(it,1,IntervalIterator);
  }
  else 
    it=it0;

#ifdef MERGE_INTERVAL_ORIENTATIONS
  if (start<0) { /* NEED TO CONVERT TO POSITIVE ORIENTATION */
    j=start;
    start= -end;
    end= -j;
    ori_sign = -1;
  }
#endif
  if (it->n == 0) { /* DEFAULT: SEARCH THE TOP NESTED LIST */
    it->n=n;
    it->i=find_overlap_start(start,end,im,n);
  }
  do {
    while (it->i>=0 && it->i<it->n && HAS_OVERLAP_POSITIVE(im[it->i],start,end)) {
      memcpy(buf+ibuf,im + it->i,sizeof(IntervalMap)); /*SAVE THIS HIT TO BUFFER */
      ibuf++;
      k=im[it->i].sublist; /* GET SUBLIST OF i IF ANY */
      it->i++; /* ADVANCE TO NEXT INTERVAL */
      if (k>=0 && (j=find_suboverlap_start(start,end,k,im,subheader,nlists))>=0) {
	PUSH_ITERATOR_STACK(it,it2,IntervalIterator); /* RECURSE TO SUBLIST */
	it2->i = j; /* START OF OVERLAPPING HITS IN THIS SUBLIST */
	it2->n = subheader[k].start+subheader[k].len; /* END OF SUBLIST */
	it=it2; /* PUSH THE ITERATOR STACK */
      }
      if (ibuf>=nbuf) /* FILLED THE BUFFER, RETURN THE RESULTS SO FAR */
	goto finally_return_result;
    }
  } while (POP_ITERATOR_STACK(it));  /* IF STACK EXHAUSTED,  EXIT */
  if (!it0) /* FREE THE ITERATOR WE CREATED.  NO NEED TO RETURN IT TO USER */
    free_interval_iterator(it);
  it=NULL;  /* ITERATOR IS EXHAUSTED */

 finally_return_result:  
#ifdef MERGE_INTERVAL_ORIENTATIONS
  reorient_intervals(ibuf,buf,ori_sign); /* REORIENT INTERVALS TO MATCH QUERY ORI */
#endif
  *p_nreturn=ibuf; /* #INTERVALS FOUND IN THIS PASS */
  return it; /* HAND BACK ITERATOR FOR CONTINUING THE SEARCH, IF ANY */
 handle_malloc_failure:
  return NULL;
}





/****************************************************************
 *
 *   FILE-BASED SEARCH FUNCTIONS
 */


/* READ A BLOCK FROM THE DATABASE FILE */
int read_imdiv(FILE *ifile,IntervalMap imdiv[],int div,int i_div,int ntop)
{
  int block;
  long ipos;
  ipos=div*i_div; /* CALCULATE POSITION IN RECORDS */
  if (ipos+div<=ntop) /* GET A WHOLE BLOCK */
    block=div;
  else /* JUST READ PARTIAL BLOCK AT END */
    block=ntop%div;
  ipos *= sizeof(IntervalMap); /* CALCULATE FILE POSITION IN BYTES */
  fseek(ifile,ipos,SEEK_SET);
  fread(imdiv,sizeof(IntervalMap),block,ifile);
  return block;
}


/* READ A SUBLIST FROM DATABASE FILE */
IntervalMap *read_sublist(FILE *ifile,SublistHeader subheader[],int isub)
{
  long ipos;
  IntervalMap *im=NULL;
  CALLOC(im,subheader[isub].len,IntervalMap);
  ipos=subheader[isub].start;
  ipos*=sizeof(IntervalMap);
  fseek(ifile,ipos,SEEK_SET);
  fread(im,sizeof(IntervalMap),subheader[isub].len,ifile);
  return im;
 handle_malloc_failure:
  return NULL;
}



int find_file_start(IntervalIterator *it,int start,int end,int isub,
		    IntervalIndex ii[],int nii,
		    SublistHeader subheader[],int nlists,
		    int ntop,int div,FILE *ifile)
{
  IntervalMap *imdiv=NULL;
  int i_div= -1,offset=0,offset_div=0;
  if (isub<0)  /* TOP-LEVEL SEARCH: USE THE INDEX */
    i_div=find_index_start(start,end,ii,nii);
  else if (subheader[isub].len>div) { /* BIG SUBLIST, SO USE THE INDEX */
    offset=subheader[isub].start;
    offset_div=offset/div;/* offset GUARANTEED TO BE MULTIPLE OF div */
    ntop=subheader[isub].len;
    nii=ntop/div; /* CALCULATE SUBLIST INDEX SIZE */
    if (ntop%div) /* ONE EXTRA ENTRY FOR PARTIAL BLOCK */
      nii++;    
    i_div=find_index_start(start,end,ii+offset_div,nii);
  }

  if (i_div>=0) { /* READ A SPECIFIC BLOCK OF SIZE div */
    if (!it->im && it->n!=div) { /* NO ALLOCATION, OR WRONG SIZE */
      FREE(it->im); /* DUMP OLD ALLOCATION, AND GET A NEW CHUNK */
      CALLOC(it->im,div,IntervalMap);
    }
    it->n=read_imdiv(ifile,it->im,div,i_div+offset_div,ntop+offset);
    it->ntop=ntop+offset; /* END OF THIS LIST IN THE BINARY FILE */
    it->nii=nii+offset_div; /* SAVE INFORMATION FOR READING SUBSEQUENT BLOCKS */
    it->i_div=i_div+offset_div; /* INDEX OF THIS BLOCK IN THE BINARY FILE */
  }
  else { /* A SMALL SUBLIST: READ THE WHOLE LIST INTO MEMORY */
    FREE(it->im);  /* DUMP ANY OLD ALLOCATION */
    it->im=read_sublist(ifile,subheader,isub);
    it->n=subheader[isub].len;
    it->nii=1;
    it->i_div=0; /* INDICATE THAT THERE ARE NO ADDITIONAL BLOCKS TO READ*/
  }

  it->i=find_overlap_start(start,end,it->im,it->n);
  return it->i;
 handle_malloc_failure:
  return -1; /* NO OVERLAP */
}


IntervalIterator *find_file_intervals(IntervalIterator *it0,int start,int end,
				      IntervalIndex ii[],int nii,
				      SublistHeader subheader[],int nlists,
				      int ntop,int div,FILE *ifile,
				      IntervalMap buf[],int nbuf,
				      int *p_nreturn)
{
  IntervalIterator *it=NULL,*it2=NULL;
  int k,ibuf=0,ori_sign=1;
  if (!it0) { /* ALLOCATE AN ITERATOR IF NOT SUPPLIED*/
    CALLOC(it,1,IntervalIterator);
  }
  else 
    it=it0;

#ifdef MERGE_INTERVAL_ORIENTATIONS
  if (start<0) { /* NEED TO CONVERT TO POSITIVE ORIENTATION */
    k=start;
    start= -end;
    end= -k;
    ori_sign = -1;
  }
#endif

  if (it->n == 0)  /* DEFAULT: SEARCH THE TOP NESTED LIST */
    find_file_start(it,start,end,-1,ii,nii,subheader,nlists,ntop,div,ifile);
  
  do { /* ITERATOR STACK LOOP */
    while (it->i_div < it->nii) { /* BLOCK ITERATION LOOP */
      while (it->i>=0 && it->i<it->n /* INDIVIDUAL INTERVAL ITERATION LOOP */
	     && HAS_OVERLAP_POSITIVE(it->im[it->i],start,end)) { /*OVERLAPS!*/
	memcpy(buf+ibuf,it->im + it->i,sizeof(IntervalMap)); /*SAVE THIS HIT */
	ibuf++;
	k=it->im[it->i].sublist; /* GET SUBLIST OF i IF ANY */
	it->i++; /* ADVANCE TO NEXT INTERVAL */
	PUSH_ITERATOR_STACK(it,it2,IntervalIterator); /* RECURSE TO SUBLIST */
	if (k>=0 && find_file_start(it2,start,end,k,ii,nii,subheader,nlists,
				    ntop,div,ifile)>=0)
	  it=it2; /* PUSH THE ITERATOR STACK */
	
	if (ibuf>=nbuf)  /* FILLED THE BUFFER, RETURN THE RESULTS SO FAR */
	  goto finally_return_result;
      }
      it->i_div++; /* TRY GOING TO NEXT BLOCK */
      if (it->i == it->n  /* USED WHOLE BLOCK, SO THERE MIGHT BE MORE */
	  && it->i_div < it->nii) { /* CONTINUE TO NEXT BLOCK */
	it->n=read_imdiv(ifile,it->im,div,it->i_div,it->ntop); /*READ NEXT BLOCK*/
	it->i=0; /* PROCESS IT FROM ITS START */
      }
    }
  } while (POP_ITERATOR_STACK(it));  /* IF STACK EXHAUSTED,  EXIT */
  if (!it0) /* FREE THE ITERATOR WE CREATED.  NO NEED TO RETURN IT TO USER */
    free_interval_iterator(it);
  it=NULL;  /* ITERATOR IS EXHAUSTED */

 finally_return_result:  
#ifdef MERGE_INTERVAL_ORIENTATIONS
  reorient_intervals(ibuf,buf,ori_sign); /* REORIENT INTERVALS TO MATCH QUERY ORI */
#endif
  *p_nreturn=ibuf; /* #INTERVALS FOUND IN THIS PASS */
  return it; /* HAND BACK ITERATOR FOR CONTINUING THE SEARCH, IF ANY */
 handle_malloc_failure:
  return NULL;
}






/* FUNCTIONS FOR READING AND WRITING OF THE BINARY DATABASE FILES */

int write_padded_binary(IntervalMap im[],int n,int div,FILE *ifile)
{
  int i,npad;
  fwrite(im,sizeof(IntervalMap),n,ifile); /* SAVE THE ACTUAL DATA */
  npad=n%div;
  if (npad) {
    npad=div-npad; /* #ITEMS NEEDED TO PAD TO AN EXACT MULTIPLE OF div */
    for (i=0;i<npad;i++) /* GUARANTEED im HAS AT LEAST ONE RECORD */
      fwrite(im,sizeof(IntervalMap),1,ifile); /*THIS IS JUST PADDING */
  }
  return n+npad; /* #RECORDS SAVED */
}


int repack_subheaders(IntervalMap im[],int n,int div,
		      SublistHeader *subheader,int nlists)
{
  int i,j,*sub_map=NULL;
  SublistHeader *sub_pack=NULL;

  CALLOC(sub_map,nlists,int);
  CALLOC(sub_pack,nlists,SublistHeader);
  for (i=j=0;i<nlists;i++) { /* PLACE SUBLISTS W/ len>div AT FRONT */
    if (subheader[i].len>div) {
      memcpy(sub_pack+j,subheader+i,sizeof(SublistHeader));
      sub_map[i]=j;
      j++;
    }
  }
  for (i=0;i<nlists;i++) { /* PLACE SUBLISTS W/ len<=div AFTERWARDS */
    if (subheader[i].len<=div) {
      memcpy(sub_pack+j,subheader+i,sizeof(SublistHeader));
      sub_map[i]=j;
      j++;
    }
  }
  for (i=0;i<n;i++) /* ADJUST im[].sublist TO THE NEW LOCATIONS */
    if (im[i].sublist>=0)
      im[i].sublist=sub_map[im[i].sublist];
  memcpy(subheader,sub_pack,nlists*sizeof(SublistHeader)); /* SAVE REORDERED LIST*/

  FREE(sub_map);
  FREE(sub_pack);
  return 0;
 handle_malloc_failure:
  return -1;
}


int write_binary_index(IntervalMap im[],int n,int div,FILE *ifile)
{
  int i,j,nsave=0;
  for (i=0;i<n;i+=div) {
#ifdef MERGE_INTERVAL_ORIENTATIONS
    if (im[i].start>=0) /* FORWARD ORI */
#endif
      fwrite(&(im[i].start),sizeof(int),1,ifile);  /*SAVE start */
#ifdef MERGE_INTERVAL_ORIENTATIONS
    else { /* REVERSE ORI */
      j= - im[i].end;
      fwrite(&j,sizeof(int),1,ifile);  /*SAVE start */
    }
#endif
    j=i+div-1;
    if (j>=n)
      j=n-1;
#ifdef MERGE_INTERVAL_ORIENTATIONS
    if (im[j].start>=0)  /* FORWARD ORI */
#endif
      fwrite(&(im[j].end),sizeof(int),1,ifile);  /*SAVE end */
#ifdef MERGE_INTERVAL_ORIENTATIONS
    else { /* REVERSE ORI */
      j= - im[j].start;
      fwrite(&j,sizeof(int),1,ifile);  /*SAVE end */
    }
#endif
    nsave++;
  }
  return nsave;
}



char *write_binary_files(IntervalMap im[],int n,int ntop,int div,
			 SublistHeader *subheader,int nlists,char filestem[])
{
  int i,npad=0,nii,*subheader_pack=NULL;
  char path[2048];
  FILE *ifile=NULL,*ifile_subheader=NULL;
  SublistHeader sh_tmp;
  static char err_msg[1024];

  if(nlists>0)
	 repack_subheaders(im,n,div,subheader,nlists); /* REPACK SMALL SUBLISTS TO END */
  sprintf(path,"%s.subhead",filestem); /* SAVE THE SUBHEADER LIST */
  ifile_subheader=fopen(path,"w");
  if (!ifile_subheader) {
    sprintf(err_msg,"unable to open file %s for writing",path);
    return err_msg;
  }
  sprintf(path,"%s.idb",filestem); /* SAVE THE DATABASE */
  ifile=fopen(path,"w");
  if (!ifile) {
    sprintf(err_msg,"unable to open file %s for writing",path);
    return err_msg;
  }
  npad=write_padded_binary(im,ntop,div,ifile); /* WRITE THE TOP LEVEL LIST */
  for (i=0;i<nlists;i++) {
    sh_tmp.start=npad; /* FILE LOCATION WHERE THIS SUBLIST STORED */
    sh_tmp.len=subheader[i].len; /* SAVE THE TRUE SUBLIST LENGTH, UNPADDED */
    fwrite(&sh_tmp,sizeof(SublistHeader),1,ifile_subheader);
    if (subheader[i].len>div) /* BIG LIST: PAD TO EXACT MULTIPLE OF div */
      npad+=write_padded_binary(im+subheader[i].start,subheader[i].len,div,ifile);
    else { /* SMALL LIST: SAVE W/O PADDING */
      fwrite(im+subheader[i].start,sizeof(IntervalMap),subheader[i].len,ifile);
      npad+=subheader[i].len;
    }
  }
  fclose(ifile);
  fclose(ifile_subheader);

  sprintf(path,"%s.index",filestem); /* SAVE THE COMPACTED INDEX */
  ifile=fopen(path,"w");
  if (!ifile) {
    sprintf(err_msg,"unable to open file %s for writing",path);
    return err_msg;
  }
  nii=write_binary_index(im,ntop,div,ifile);
  for (i=0;i<nlists;i++) /* ALSO STORE INDEX DATA FOR BIG SUBLISTS */
    if (subheader[i].len>div)
      nii+=write_binary_index(im+subheader[i].start,subheader[i].len,div,ifile);
  fclose(ifile);

  sprintf(path,"%s.size",filestem); /* SAVE BASIC SIZE INFO*/
  ifile=fopen(path,"w");
  if (!ifile) {
    sprintf(err_msg,"unable to open file %s for writing",path);
    return err_msg;
  }
  fprintf(ifile,"%d %d %d %d %d\n",n,ntop,div,nlists,nii);
  fclose(ifile);

  return NULL; /* RETURN CODE SIGNALS SUCCESS!! */
}



IntervalDBFile *read_binary_files(char filestem[],char err_msg[])
{
  int n,ntop,div,nlists,nii;
  char path[2048];
  IntervalIndex *ii=NULL;
  SublistHeader *subheader=NULL;
  IntervalDBFile *idb_file=NULL;
  FILE *ifile=NULL;

  sprintf(path,"%s.size",filestem); /* SAVE BASIC SIZE INFO*/
  ifile=fopen(path,"r");
  if (!ifile) {
    if (err_msg)
      sprintf(err_msg,"unable to open file %s",path);
    return NULL;
  }
  fscanf(ifile,"%d %d %d %d %d",&n,&ntop,&div,&nlists,&nii);
  fclose(ifile);

  CALLOC(ii,nii,IntervalIndex);
  sprintf(path,"%s.index",filestem); /* SAVE BASIC SIZE INFO*/
  ifile=fopen(path,"r");
  if (!ifile) {
    if (err_msg)
      sprintf(err_msg,"unable to open file %s",path);
    return NULL;
  }
  fread(ii,sizeof(IntervalIndex),nii,ifile);
  fclose(ifile);

  if(nlists>0){
	 CALLOC(subheader,nlists,SublistHeader);
	 sprintf(path,"%s.subhead",filestem); /* SAVE THE SUBHEADER LIST */
	 ifile=fopen(path,"r");
	 if (!ifile) {
		if (err_msg)
		  sprintf(err_msg,"unable to open file %s",path);
		return NULL;
	 }
	 fread(subheader,sizeof(SublistHeader),nlists,ifile);  /*SAVE LIST */
	 fclose(ifile);
  }

  CALLOC(idb_file,1,IntervalDBFile);
  idb_file->n=n;
  idb_file->ntop=ntop;
  idb_file->nlists=nlists;
  idb_file->div=div;
  idb_file->nii=ntop/div;
  if (ntop%div) /* INDEX IS PADDED TO EXACT MULTIPLE OF div */
    idb_file->nii++; /* ONE EXTRA ENTRY FOR PARTIAL BLOCK */
  idb_file->ii=ii;
  idb_file->subheader=subheader;
  sprintf(path,"%s.idb",filestem); /* SAVE BASIC SIZE INFO*/
  idb_file->ifile_idb=fopen(path,"r");
  if (!idb_file->ifile_idb) {
    if (err_msg)
      sprintf(err_msg,"unable to open file %s",path);
    free(idb_file);
    return NULL;
  }
  return idb_file;
 handle_malloc_failure:
  FREE(ii); /* DUMP OUR MEMORY */
  FREE(subheader);
  FREE(idb_file);
  return NULL;
}



int free_interval_dbfile(IntervalDBFile *db_file)
{
  fclose(db_file->ifile_idb);
  FREE(db_file->ii);
  FREE(db_file->subheader);
  free(db_file);
  return 0;
}
