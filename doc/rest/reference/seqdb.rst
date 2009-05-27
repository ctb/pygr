:mod:`seqdb` --- Sequence database interfaces
=============================================

.. module:: seqdb
   :synopsis: Sequence database interfaces.
.. moduleauthor:: Christopher Lee <leec@chem.ucla.edu>
.. sectionauthor:: Christopher Lee <leec@chem.ucla.edu>


The seqdb module provides a simple, consistent interface to sequence databases 
from a variety of different storage sources such as FASTA, BLAST and 
relational databases.  Sequence databases are modeled 
(like other Pygr container classes) as dictionaries, whose keys are 
sequence IDs and whose values are sequence objects.  
Pygr sequence objects use the Python sequence protocol in all the 
ways you'd expect: a subinterval of a sequence object is just a 
Python slice (s[0:10]), which just returns a sequence object 
representing that interval; the reverse complement is just -s; 
the length of a sequence is just len(s); to obtain the actual 
string sequence of a sequence object is just str(s).  
Pygr sequence objects work intelligently with different 
types of back-end storage (e.g. relational databases or BLAST databases) 
to efficiently access just the parts of sequence that are requested, 
only when an actual sequence string is needed.

SequenceFileDB
--------------
Interface to a sequence database, typically initialized from a FASTA sequence file.
This class was renamed from "BlastDB" in Pygr 0.8.

Options for constructing a SequenceFileDB:

.. class:: SequenceFileDB(filepath=None, itemClass=FileDBSequence, itemSliceClass=None, reader=None, autoGC=True, **kwargs)

   Open a sequence file as a "database" object, giving the user access
   to its sequences.

   * *filepath*: path to the text sequence file (typically FASTA).

   * *itemClass*: the object class to use for instantiating new sequence
     objects from this database.  You can set this to create customized
     sequence behaviors.
     This is used by :class:`pygr.Data` to propagate correct attribute schemas to
     items / slices from database containers managed by it.

   * *itemSliceClass*: the object class to use for instantiating new
     sequence slice objects (i.e. subintervals of sequences from this database).
     You can set this to create customized sequence behaviors.
     This is used by :class:`pygr.Data` to propagate correct attribute schemas to
     items / slices from database containers managed by it.

   * *reader*: allows you to specify a parser function.
     It will be called with
     two arguments: ``reader(ifile, filename)``; and it should
     act like a generator that yields one or more objects that
     each provide a sequence ID, length and sequence string.  See
     the Pygr Developer Guide for details.

   * *autoGC*: if True, automatically garbage-collect cache entries
     for sequence objects.  When you request a sequence object from
     SequenceFileDB, it keeps it in a cache, both for speeding up future
     requests, and to ensure that any two requests for the same sequence ID
     are guaranteed to return the same Python object.  To save memory,
     when you drop all references to a given sequence object, SequenceFileDB
     also flushes it from its cache.  This is implemented using a Python
     WeakValueDictionary.


Useful methods:

:class:`SequenceFileDB` implements a Python dictionary interface,
so all of the methods you would expect for a dictionary are available.
For example:

.. method:: iter()

   iterate over all IDs in the database.


.. method:: len()

   returns number of sequences in the database.


.. method:: __invert__()

   The invert operator (\textasciitilde, the "tilde" character)
   enables reverse-mapping of sequence objects to their string ID::

      id = (~db)[seq] # GET IDENTIFIER FOR THIS SEQUENCE FROM ITS DATABASE





Useful attributes:

  
.. attribute:: filepath

   the location of the raw sequence file (by default, FASTA)
   upon which this :class:`SequenceFileDB` is based.
  


PrefixUnionDict
---------------
This class acts as a wrapper for a set of dictionaries, each
of which is assigned a specific string "prefix".  It provides
a dictionary interface that accepts string keys of the form
"prefix.suffix", and returns d['suffix'] where *d* is
the dictionary associated with the corresponding prefix.  This
is useful for providing a unified "database interface" to a
set of multiple databases::

   hg17 = BlastDB('/usr/tmp/ucsc_msa/hg17')
   mm5 = BlastDB('/usr/tmp/ucsc_msa/mm5')

   genomeUnion = PrefixUnionDict({ 'mm5': mm5, 'hg17', hg17 })

   # retrieve a mouse chr
   mouse_chr3 = genomeUnion['mm5.chr3']
   
   # retrieve a human chr
   human_chr2 = genomeUnion['hg17.chr2']

   # check to see if sequence ival is in genomeUnion:
   if human_chr2[5000:6000] in genomeUnion: ...

It provides a :meth:`__contains__` method that tests whether
a given sequence object is derived from the :class:`PrefixUnionDict`
(see example above).  Here are some additional methods:

.. class:: PrefixUnionDict(prefixDict=None, separator='.', filename=None, dbClass=SequenceFileDB, trypath=None)

   You can create a :class:`PrefixUnionDict` either using
   a *prefixDict* (whose keys are string prefixes, and whose
   values are sequence databases), or using a previously created
   header file *filename*.
   Using the header file, the constructor will
   automatically open all the sequence databases for you.
   When opening from a header file, you can also specify a
   *dbClass* to be used for opening individual sequence databases
   listed in the header file; the default is :class:`SequenceFileDB`.
   The database class constructor must take a single argument,
   which is the filepath for opening the database.  The
   *separator* character is used to form "prefix.suffix"
   identifiers.  *trypath* is a list of paths to search for
   sequence dbs to open.


.. method:: __invert__()

   The invert operator (~, the "tilde" character) enables
   reverse-mapping of sequence objects to their string ID, that is,
   for a given sequence object *seq* derived from the union (or a
   slice of a sequence from the union), return a string identifier in
   the form of "foo.bar".

   This is the recommended way to get the "fully qualified sequence
   ID", i.e. with the appropriate prefix prepended::

      id = (~db)[seq]
      # e.g. 'hg17.chr3'

.. method:: newMemberDict(**kwargs)

   Returns a new member dictionary for testing membership in
   the distinct prefix groups.  See :class:`PrefixUnionMemberDict`.


.. method:: cacheHint(owner,ivalDict)

   Communicates a set of caching hints to the appropriate member
   databases.  *ivalDict* must be a dictionary whose keys are
   sequence ID strings, and whose values are each a (start,stop) tuple
   for the associated covering interval coordinates to
   cache for each sequence.  *owner* should be a python object
   whose existence controls the lifetime of these cache hints.
   When *owner* is garbage-collected by Python (due to its
   reference count going to zero), the member databases will clear
   these cache hints from their cache storage.

   On :class:`PrefixUnionDict`, this method simply passes along
   the cache hints to the appropriate member databases by calling
   their :meth:`cacheHint` method, without itself doing anything
   to cache the information.

.. method:: getName(path)

   **This method is deprecated;** instead use the :meth:`__invert__` operator
   above.


.. method:: writeHeaderFile(filename)

   **This method is deprecated,** because it is restricted to assuming
   that all sequence dictionaries it contains are of a single class.
   We recommend that you instead save it to pygr.Data, or pickle it
   directly using pygr.Data.dumps().

   Save a header file for this union, to reopen later.  It saves the
   separator character, and a list of prefixes and filepaths to the
   various sequence databases (which must have a :attr:`filepath`
   attribute).  This header file can be used for later reopening the
   prefix-union in a single step.


_PrefixUnionMemberDict
----------------------
Implements membership testing on distinct prefix groups.  Specifically,
you can bind a given prefix to a value::

   d['prefix1'] = value

then test whether a given object *k* is a member of any of the
prefix groups in the dictionary::

   v = d[k] # raises KeyError if k not a member of 'prefix1' or other prefix group in d


.. class:: _PrefixUnionMemberDict(puDict,default=None,attrMethod=lambda x:x.pathForward.db)

   * *puDict* must be a :class:`PrefixUnionDict`, whose prefix groups constitute the
     allowed possible key groups for this membership dictionary.  *default*
     provides a default value to apply to any key whose prefix has not been explicitly
     given a value in this dictionary.  If no *default* is set, this dictionary
     will raise a :exc:`KeyError` for any key whose prefix has not been
     explicitly given a value in this dictionary.

   * *attrMethod* specifies a method for obtaining
     the actual prefix container object from a given member key object.  The default

   * *attrMethod* treats the key as a sequence object and tries to determine what
     database container it is from.


.. method:: possibleKeys()

   Returns an iterator for the key values (prefix strings) that are allowed for
   this dictionary, obtained from the bound :class:`PrefixUnionDict`.


PrefixDictInverse
-----------------
Provides the interface to the inverse mapping of the :class:`PrefixUnionDict`.

.. method:: __getitem__(k)

   Returns the fully-qualified string ID for sequence object *k*.
   Properly handles both sequence annotation object and regular sequence
   objects.


PrefixDictInverseAdder
----------------------
Adds the capability of automatically adding new sequence databases to the
:class:`PrefixUnionDict`, if needed.  This is implemented by extending
the standard :meth:`__getitem__` method:
.. method:: __getitem__(k)

   Returns the fully-qualified string ID for sequence object *k*.
   Properly handles both sequence annotation object and regular sequence
   objects.  If sequence object *k* is from a sequence database that
   is not in the :class:`PrefixUnionDict`, it will be automatically added
   to the prefixUnion, if the prefixUnion has an :attr:`addAll` attribute
   set to *True*; if not, a :exc:`KeyError` is raised.
   This is used in the standard :class:`NLMSA` write mode 'w'
   to allow users to add sequences to the alignment without having to
   previously add the sequence databases containing those sequences,
   to the prefixUnion for the NLMSA.




FileDBSequence
--------------
The default class for sequence objects returned from SequenceFileDB.
It provides efficient, fast access to sequence slices (subsequences).
When the SequenceFileDB is initially opened,
Pygr constructs a length and offset index that enables FileDBSequence to ``seek()``
to the correct location for any substring of the sequence.

SQLSequence
-----------

Implements a subclass inheriting from SQLRow and SequenceBase, to use a relational database table to obtain the actual sequence.  There are three minor variants DNASQLSequence, RNASQLSequence, ProteinSQLSequence (so that the sequence does not have to analyze itself to determine what kind of sequence it is).  Its constructor takes the same arguments as SQLRow(table, id), where table is the SQLTable object representing the table in which the sequence is stored, and id is the primary key of the row representing this sequence.  However, normally this class is simply passed to the Table object itself so that it will use it to instantiate new row objects whenever they are requested via its dictionary interface.

*Python DB-API 2.0*: this class conforms to the Python DB-API 2.0.
Typically you must supply a DB-API 2.0-compliant database cursor to the
:class:`SQLTable` constructor.  To do so, you must have some DB-API 2.0-compliant
module (such as :mod:`MySQLdb`) installed for connecting to a database server.

Here's a simple example of customizing SQLSequence for your data::

   class YiProteinSequence(ProteinSQLSequence): # CREATE A NEW SQL SEQUENCE CLASS
       def __len__(self): return self.protein_length  # USE LENGTH STORED IN DATABASE
   protein = jun03[protein_seq_t] # protein IS OUR SQLTable OBJECT REPRESENTING PROTEIN SEQUENCE TABLE
   protein.objclass(YiProteinSequence) # FORCE PROTEIN SEQ TABLE TO USE THIS TO INSTANTIATE ROW OBJECTS
   pseq = protein['Hs.1162'] # GET PROTEIN SEQUENCE OBJECT FOR A SPECIFIC CLUSTER


Let's go through this line by line:


  
* we create a subclass of ProteinSQLSequence to show how Python makes it easy to create customized behaviors that can make database access more efficient.  Here we've simply added a __len__ method that uses the protein_length attribute obtained directly from the database, courtesy of SQLRow.__getattr__, which knows what columns exist in the database, and provides them transparently as object attributes.  (The ordinary SequenceBase __len__ method calculates it by obtaining the whole sequence string and calculating its length.  Clearly it's more efficient for the database to retrieve this number (stored as a column called protein_length) and return it, rather than making it send us the whole sequence).
  
* next we call the protein.objclass() method to inform the table object that it should use our new class for instantiating any row objects for this table.  It will call this class with the usual SQLRow contructor arguments (table, id).


BlastDBXMLRPC
-------------
A subclass of :class:`SequenceFileDB` that adds a couple methods needed to serve
the data to clients connecting over XMLRPC.  For example, to make an XMLRPC
server for a blast database, accessible on port 5020::

   import coordinator
   server = coordinator.XMLRPCServerBase(name,port=5020)
   db = BlastDBXMLRPC('sp') # OPEN BlastDB AS USUAL, BUT WITH SUBCLASS
   server['sp'] = db # ADD OUR DATABASE TO THE XMLRPC SERVER
   server.serve_forever() # START SERVING XMLRPC REQUESTS, UNTIL KILLED.


XMLRPCSequenceDB
----------------
Class for a client interface that accesses a Blast database over
XMLRPC (from the the :class:`BlastDBXMLRPC` acting as the server).

.. class:: XMLRPCSequenceDB(url,name)

   *url* must be the URL (including port number) for accessing the
   XMLRPC server; *name* must be the key of the BlastDBXMLRPC object
   in that server's dictionary (in the example above, it would be 'sp').
   Thus to access the server above (assuming it is running on leelab port 5020)::

      db = XMLRPCSequenceDB('http://leelab:5020','sp')
      hbb = db['HBB_HUMAN'] # GET A SEQUENCE OBJECT FROM THE DATABASE...


Currently, this class provides sequence access.  You can work with sequences
exactly as you would with a :class:`BlastDB`, but cannot perform actual BLAST searches
(i.e. the :meth:`blast` and :meth:`megablast` methods don't work over XMLRPC).
