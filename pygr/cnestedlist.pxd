
cdef extern from "string.h":
  ctypedef int size_t
  void *memcpy(void *dst,void *src,size_t len)
  void *memmove(void *dst,void *src,size_t len)

cdef extern from "stdlib.h":
  void free(void *)
  void *realloc(void *,size_t)
  int c_abs "abs" (int)
  void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))

cdef extern from "stdio.h":
  ctypedef struct FILE:
    pass
  FILE *fopen(char *,char *)
  int fclose(FILE *)
  int sscanf(char *str,char *fmt,...)
  int sprintf(char *str,char *fmt,...)


cdef extern from "intervaldb.h":
  ctypedef struct IntervalMap:
    int start
    int end
    int target_id
    int target_start
    int target_end
    int sublist

  ctypedef struct IntervalIndex:
    int start
    int end

  ctypedef struct SublistHeader:
    int start
    int len
  
  ctypedef struct IntervalDBFile:
    int n
    int ntop
    int nlists
    int div
    int nii
    IntervalIndex *ii
    SublistHeader *subheader
    FILE *ifile_idb

  ctypedef struct IntervalIterator:
    pass

  ctypedef struct IDInterval:
    int id
    int start
    int stop

  int imstart_qsort_cmp(void *void_a,void *void_b)
  IntervalMap *read_intervals(int n,FILE *ifile)
  SublistHeader *build_nested_list(IntervalMap im[],int n,int *p_n,int *p_nlists)
  IntervalMap *interval_map_alloc(int n)
  IntervalIterator *interval_iterator_alloc()
  int free_interval_iterator(IntervalIterator *it)
  IntervalIterator *reset_interval_iterator(IntervalIterator *it)
  IntervalIterator *find_intervals(IntervalIterator *it0,int start,int end,IntervalMap im[],int n,SublistHeader subheader[],int nlists,IntervalMap buf[],int nbuf,int *p_nreturn)
  char *write_binary_files(IntervalMap im[],int n,int ntop,int div,SublistHeader *subheader,int nlists,char filestem[])
  IntervalDBFile *read_binary_files(char filestem[],char err_msg[])
  int free_interval_dbfile(IntervalDBFile *db_file)
  IntervalIterator *find_file_intervals(IntervalIterator *it0,int start,int end,IntervalIndex ii[],int nii,SublistHeader subheader[],int nlists,int ntop,int div,FILE *ifile,IntervalMap buf[],int nbuf,int *p_nreturn)
  int write_padded_binary(IntervalMap im[],int n,int div,FILE *ifile)
  int read_imdiv(FILE *ifile,IntervalMap imdiv[],int div,int i_div,int ntop)
  IDInterval *interval_id_alloc(int n)
  int interval_id_union(int id,int start,int stop,IDInterval iv[],int n)
  IDInterval *interval_id_compact(IDInterval iv[],int *p_n)



cdef class IntervalDB:
  cdef int n
  cdef int ntop
  cdef int nlists
  cdef IntervalMap *im
  cdef SublistHeader *subheader


cdef class IntervalDBIterator:
  cdef IntervalIterator *it,*it_alloc
  cdef IntervalMap im_buf[1024]
  cdef int ihit,nhit,start,end
  cdef IntervalDB db

  cdef int cnext(self)

cdef class IntervalFileDB:
  cdef IntervalDBFile *db

cdef class IntervalFileDBIterator:
  cdef IntervalIterator *it,*it_alloc
  cdef IntervalMap *im_buf
  cdef int ihit,nhit,start,end,nbuf
  cdef IntervalFileDB db

  cdef int cnext(self,int *pkeep)
  cdef int extend(self,int ikeep)
  cdef int saveInterval(self,int start,int end,int target_id,
                        int target_start,int target_end)
  cdef int nextBlock(self,int *pkeep)
  cdef IntervalMap *getIntervalMap(self)
  cdef int loadAll(self)


cdef class NLMSALetters:
  cdef readonly object pathstem
  cdef readonly object seqs
  cdef readonly object seqlist
  cdef int do_build
  cdef int lpo_id

  cdef int is_lpo(self,int id)

cdef class NLMSASequence:
  cdef readonly int id,length,nbuild,fixed_length
  cdef readonly object seq
  cdef IntervalFileDB db
  cdef FILE *build_ifile
  cdef readonly object filestem
  cdef NLMSALetters nlmsaLetters

cdef class NLMSASlice:
  cdef readonly start,stop 
  cdef int n,nseqBounds
  cdef IntervalMap *im
  cdef IDInterval *seqBounds
  cdef NLMSASequence nlmsaSequence


cdef class NLMSANode:
  cdef readonly int id,ipos
  cdef int istart,istop,n
  cdef NLMSASlice nlmsaSlice

  cdef int check_edge(self,int iseq,int ipos)


cdef class NLMSASliceIterator:
  cdef int ipos,istart,istop
  cdef NLMSASlice nlmsaSlice

  cdef int advanceStartStop(self)