import unittest
from testlib import testutil, PygrTestProgram
from pygr import sequence, seqdb

class Sequence_Test(unittest.TestCase):
    'basic sequence class tests'
    
    def setUp(self):
        self.seq = sequence.Sequence('atttgactatgctccag', 'foo')     
    
    def test_length(self):
        "Sequence lenght"
        assert len(self.seq) == 17
    
    def test_slice(self):
        "Sequence slice"
        assert str(self.seq[5:10]) == 'actat'
    
    def test_slicerc(self):
        "Sequence slice then reverse complement"
        assert str(-(self.seq[5:10])) == 'atagt'
    
    def test_rcslice(self):
        "Sequence reverse complement then slice"
        assert str((-self.seq)[5:10]) == 'gcata'
    
    def test_truncate(self):
        "Sequence truncate"
        assert str(self.seq[-202020202:5]) == 'atttg'
        assert self.seq[-202020202:5] == self.seq[0:5]
        assert self.seq[-2020202:] == self.seq
        assert str(self.seq[-202020202:-5]) == 'atttgactatgc'
        assert str(self.seq[-5:2029]) == 'tccag'
        assert str(self.seq[-5:]) == 'tccag'
        try:
            self.seq[999:10000]
            raise ValueError('failed to trap out of bounds slice')
        except IndexError:
            pass
        try:
            self.seq[-10000:-3000]
            raise ValueError('failed to trap out of bounds slice')
        except IndexError:
            pass
        try:
            self.seq[1000:]
            raise ValueError('failed to trap out of bounds slice')
        except IndexError:
            pass
    
    def test_rctruncate(self):
        "Sequence reverse complement truncate"
        seq= -self.seq
        assert str(seq[-202020202:5]) == 'ctgga'
        assert seq[-202020202:5] == seq[0:5]
        assert seq[-2020202:] == seq
        assert str(seq[-202020202:-5]) == 'ctggagcatagt'
        assert str(seq[-5:2029]) == 'caaat'
        assert str(seq[-5:]) == 'caaat'
        try:
            seq[999:10000]
            raise ValueError('failed to trap out of bounds slice')
        except IndexError:
            pass
        try:
            seq[-10000:-3000]
            raise ValueError('failed to trap out of bounds slice')
        except IndexError:
            pass
        try:
            seq[1000:]
            raise ValueError('failed to trap out of bounds slice')
        except IndexError:
            pass
    
    def test_join(self):
        "Sequence join"
        assert str(self.seq[5:15] * self.seq[8:]) == 'atgctcc'

    def test_rcjoin(self):
        "Sequence reverse complement join"
        assert str((-(self.seq[5:10])) * ((-self.seq)[5:10])) == 'ata'

    def test_seqtype(self):
        "Sequence lenght"
        assert self.seq.seqtype() == sequence.DNA_SEQTYPE
        assert sequence.Sequence('auuugacuaugcuccag', 'foo').seqtype() == \
                         sequence.RNA_SEQTYPE
        assert sequence.Sequence('kqwestvvarphal', 'foo').seqtype() == \
                         sequence.PROTEIN_SEQTYPE

class SequenceTranslation_Test(unittest.TestCase):
    '''Tests for the translation() function.'''
    
    def setUp(self):
        self.db = seqdb.SequenceFileDB(testutil.datafile('translation.fa'))
        self.seq6 = self.db['seq6']
        self.seq7 = self.db['seq7']
        self.seq8 = self.db['seq8']
        self.seq9 = self.db['seq9']

    def test_6_path(self):
        assert len(self.seq6.translation(1)) == 2
        assert len(self.seq6.translation(2)) == 1
        assert len(self.seq6.translation(3)) == 1
        
    def test_7_path(self):
        assert len(self.seq7.translation(1)) == 2
        assert len(self.seq7.translation(2)) == 2
        assert len(self.seq7.translation(3)) == 1
        
    def test_8_path(self):
        assert len(self.seq8.translation(1)) == 2
        assert len(self.seq8.translation(2)) == 2
        assert len(self.seq8.translation(3)) == 2

    def test_slice_frame_adjust(self):
        s = self.seq8
        assert str(s[1:].translation(1)) == str(s[0:].translation(2))
        assert str(s[2:].translation(1)) == str(s[0:].translation(3))
        assert str(s[3:].translation(1)) == str(s[0:].translation(1))

    def test_slice_frame_adjust(self):
        s = self.seq9

        ### slices of first-frame translation are correct:
        
        assert str(s[0:3].translation(1)) == 'M', s[0:3].translation(1)
        assert str(s[3:6].translation(1)) == '*', s[3:6].translation(1)
        assert str(s[6:9].translation(1)) == 'L', s[6:9].translation(1)

        ### starting with nt 0
        
        x = s[0:3].translation(1)                    # frame 1 @ nt 0
        y = s.translation(1)[:1]                     # == frame 1, aa 0
        assert str(x) == str(y), (str(x), str(y))
        assert x.frame == 1
        assert y.frame == 1

        x = s[0:4].translation(2)                    # frame 2 @ nt 0
        y = s.translation(2)[:1]                     # == frame 2, aa 0
        assert str(x) == str(y), (str(x), str(y))
        assert x.frame == 2
        assert y.frame == 2

        x = s[0:5].translation(3)                    # frame 3 @ nt 0
        y = s.translation(3)[:1]                     # == frame 3, aa 0
        assert str(x) == str(y), (str(x), str(y))
        assert x.frame == 3
        assert y.frame == 3

        ## starting with nt 1

        x = s[1:4].translation(1)                    # frame 1 @ nt 1
        y = s.translation(2)[:1]                     # == frame 2, aa 0
        assert str(x) == str(y)
        assert y.frame == 2

        x = s[1:5].translation(2)                    # frame 2 @ nt 1
        y = s.translation(3)[:1]                     # == frame 3, aa 0
        assert str(x) == str(y), (str(x), str(y))
        assert y.frame == 3

        x = s[1:6].translation(3)                    # frame 3 @ pos 1
        y = s.translation(1)[1:2]                    # == frame 1, aa 1
        assert str(x) == str(y), (str(x), str(y))
        assert y.frame == 1

        ## starting with nt 2

        x = s[2:5].translation(1)                    # frame 1 @ nt2
        y = s.translation(3)[0]                      # == frame 3, aa 0
        assert str(x) == str(y), (str(x), str(y))
        assert y.frame == 3

        x = s[2:6].translation(2)                    # frame 2 @ nt2
        y = s.translation(1)[1:2]                    # == frame 1, aa 1
        assert str(x) == str(y)
        assert y.frame == 1

        x = s[2:7].translation(3)                    # frame 3 @ nt2
        y = s.translation(2)[1:2]                    # == frame 2, aa 1
        assert str(x) == str(y)
        assert y.frame == 2

    def test_annot_identity(self):
        return
    
        s = self.seq8
        
        x = s[1:3].translation(1)
        y = s.translation(2)[:1]
        assert repr(x) == repr(y), (x, y)           # #@CTB

        print '----'

        x = s[1:4].translation(2)
        y = s.translation(3)[:1]
        
        print x, x.frame, y, y.frame
        print '----'
        assert repr(x) == repr(y), (x, y)           # #@CTB

# @CTB
'''
#from pygrdata_test import PygrSwissprotBase
class Blast_Test(PygrSwissprotBase):
    'test basic blast functionality'
    @skip_errors(OSError,KeyError)
    def setup(self):
        PygrSwissprotBase.setup(self)
        import pygr.Data
        self.sp = pygr.Data.Bio.Seq.Swissprot.sp42()
        import os
        blastIndexPath = os.path.join(os.path.dirname(self.sp.filepath),
                                      'wikiwacky')
        self.sp.formatdb(blastIndexPath)
    def blast(self):
        hbb = self.sp['HBB1_TORMA']
        hits = self.sp.blast(hbb)
        edges = hits[hbb].edges(maxgap=1,maxinsert=1,
                                minAlignSize=14,pIdentityMin=0.5)
        for t in edges:
            assert len(t[0])>=14, 'result shorter than minAlignSize!'
        result = [(t[0],t[1],t[2].pIdentity()) for t in edges]
        store = PygrDataTextFile(os.path.join('results', 'seqdb1.pickle'))
        correct = store['hbb blast 1']
        assert approximate_cmp(result,correct,.0001) == 0, 'blast results should match'
        result = [(t[0],t[1],t[2].pIdentity()) for t in hits[hbb].generateSeqEnds()]
        correct = store['hbb blast 2']
        assert approximate_cmp(result,correct,.0001) == 0, 'blast results should match'
        trypsin = self.sp['PRCA_ANASP']
        try:
            hits[trypsin]
            raise ValueError('failed to catch bad alignment query')
        except KeyError:
            pass
class Blast_reindex_untest(Blast_Test):
    'test building blast indexes under a different name'
    @skip_errors(OSError,KeyError)
    def setup(self):
        PygrSwissprotBase.setup(self)
        import pygr.Data
        self.sp = pygr.Data.Bio.Seq.Swissprot.sp42()
        import os
        blastIndexPath = os.path.join(os.path.dirname(self.sp.filepath),'wikiwacky')
        self.sp.formatdb()
        #self.sp.formatdb(blastIndexPath) # FORCE IT TO STORE INDEX WITH DIFFERENT NAME
        #print 'blastIndexPath is',self.sp.blastIndexPath

'''

if __name__ == '__main__':
    PygrTestProgram(verbosity=2)
