Creating, Querying, and Storing Alignments
------------------------------------------

Alignment Basics
^^^^^^^^^^^^^^^^

Pygr multiple alignment objects can be treated as mappings of sequence
intervals onto sequence intervals.  Here is an example showing basic
operations for constructing and querying an alignment.

First, create an empty in-memory alignment:

   >>> from pygr import cnestedlist
   >>> simple_al = cnestedlist.NLMSA('hbb', mode='memory')

Load some sequences, too:

   >>> from pygr import seqdb
   >>> db = seqdb.SequenceFileDB('data/sp_all_hbb')
   >>> mouse = db['HBB1_MOUSE']
   >>> rat = db['HBB1_RAT']
   >>> frog = db['HBB1_XENLA']

Now, add the mouse sequence to the alignment:

   >>> simple_al += mouse

and align several intervals from other sequences to the mouse sequence:

   >>> ival = mouse[40:60]
   >>> simple_al[ival] += rat[42:62]
   >>> simple_al[ival] += frog[38:58]

Once we're done adding aligned intervals, build the alignment object
to prepare it for querying:

   >>> simple_al.build()

Now we can query the alignment with any new interval overlapping the
mouse interval to which we aligned everything:

   >>> sub_ival = mouse[48:52]
   >>> for aligned_ival in simple_al[sub_ival]:
   ...   print repr(aligned_ival)
   HBB1_RAT[50:54]
   HBB1_XENLA[46:50]

We can also use any of the other sequences as keys, e.g. the whole
rat sequence:

   >>> for aligned_ival in simple_al[rat]:
   ...   print repr(aligned_ival)
   HBB1_MOUSE[40:60]

When you query 'simple_al' with an interval, it returns an
``NLMSASlice`` object.  The examples above all use the basic iterator
interface to ``NLMSASlice``, and it's equivalent to calling ``keys()``
with no arguments (we'll discuss ``keys()`` in more detail below).
You can also call ``edges()``, which will return a triple of source
interval, destination interval, and edge information:

   >>> for src, dest, edge in simple_al[mouse].edges():
   ...   print repr(src), 'aligns to', repr(dest)
   ...   print 'Identity across alignment:', edge.pIdentity()
   ...   print '--'
   HBB1_MOUSE[40:60] aligns to HBB1_RAT[42:62]
   Identity across alignment: 0.15
   --
   HBB1_MOUSE[40:60] aligns to HBB1_XENLA[38:58]
   Identity across alignment: 0.05
   --

You can also retrieve edges directly by querying the ``NLMSASlice`` with
the aligned sequence, treating it as a mapping:

   >>> edge = simple_al[mouse][rat]
   >>> print '%.4f' % (edge.pIdentity(),)
   0.0068

Gaps
^^^^

Storing alignments on disk
^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating an NLMSA object can take a long time and a lot of memory;
what if you want to build it just once, and then query it multiple
times?  You can do this by creating an NLMSA in 'w' (write) mode,
rather than 'memory' mode; otherwise the semantics are the same.

Create the NLMSA,

   >>> simple_al = cnestedlist.NLMSA('tempdir/hbb', mode='w')

load the sequences,

   >>> db = seqdb.SequenceFileDB('data/sp_all_hbb')
   >>> mouse = db['HBB1_MOUSE']
   >>> rat = db['HBB1_RAT']
   >>> frog = db['HBB1_XENLA']

add the mouse sequence into the alignment

   >>> simple_al += mouse

and align several intervals from other sequences to the mouse sequence:

   >>> ival = mouse[40:60]
   >>> simple_al[ival] += rat[42:62]
   >>> simple_al[ival] += frog[38:58]

And, finally, build it and delete the in-memory handle:

   >>> simple_al.build(saveSeqDict=True)
   >>> del simple_al

Now, to load this alignment, we need to specify the sequence source or
sources that we used to build it -- we can do that by using
``PrefixUnionDict`` to construct a ``seqDict`` and pass it into the NLMSA.

   >>> seqDict = seqdb.PrefixUnionDict({ 'sp_all_hbb': db })
   >>> loaded_al = cnestedlist.NLMSA('tempdir/hbb', seqDict=seqDict)
   >>> loaded_al[ival].keys()
   [HBB1_RAT[42:62], HBB1_XENLA[38:58]]

Here we can use our interval from above because the sequence references
stored in the NLMSA point to ``db``, the database that our interval came
from in the first place.

You can also load the saved seqDict (see ``simple_al.build``, above, where
we told pygr to save the sequence dictionary):

   >>> loaded_al = cnestedlist.NLMSA('tempdir/hbb')

Now, however, you can't query with our original ival, because we
loaded a new seqDict into memory. Even though it's pointing at the
same on-disk file as ``db`` did before, pygr only keeps track of
sequence-to-database relationships in memory.  So now you have to
manually retrieve the mouse sequence from the new seqDict in order to
query the NLMSA:

   >>> seqDict = loaded_al.seqDict
   >>> ival = seqDict['sp_all_hbb.HBB1_MOUSE']

and voila, now we can query the alignment, etc.

   >>> loaded_al[ival].keys()
   [HBB1_RAT[42:62], HBB1_XENLA[38:58]]

In practice, if you store your sequence collections in ``pygr.Data``,
you don't need to worry about seqDict mechanisms.  However, if you're
not using ``pygr.Data`` then you'll need to keep track of your sequence
dictionaries.

Creating alignments with BLAST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Using the "translated BLASTs" (blastx, tblastn, tblastx)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Building an Alignment Database from MAF files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example: Mapping an entire gene set onto a new genome version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
