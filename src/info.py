import _omnical as _O
import numpy as np, numpy.linalg as la
import warnings, os, time
import array # XXX need to remove this dependency
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy.sparse as sps

# XXX all this meta stuff about "info" almost assuredly means info needs to be a class
infokeys = ['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B','At','Bt','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
infokeys_optional = ['totalVisibilityId']
binaryinfokeys=['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','At','Bt']
cal_name = {0: "Lincal", 1: "Logcal"}

int_infokeys = ['nAntenna','nUBL','nBaseline']
intarray_infokeys = ['subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','A','B','At','Bt']
intarray_infokeys_optional = ['totalVisibilityId']

float_infokeys = ['antloc','ubl','degenM','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']

MARKER = 9999999

#  ___        _              _          _   ___       __     
# | _ \___ __| |_  _ _ _  __| |__ _ _ _| |_|_ _|_ _  / _|___ 
# |   / -_) _` | || | ' \/ _` / _` | ' \  _|| || ' \|  _/ _ \
# |_|_\___\__,_|\_,_|_||_\__,_\__,_|_||_\__|___|_||_|_| \___/

class RedundantInfo(_O.RedundantInfo):
    '''Container for metadata used by redundant calibrator, which is eventually passed into C++ routines.'''
    def __init__(self, filename=None, verbose=False, preview_only=False, txtmode=False, threshold=128):
        _O.RedundantInfo.__init__(self)
        self.threshold = threshold # XXX move to versions
        if filename:
            if txtmode: self.fromfile_txt(filename, verbose=verbose, preview_only=preview_only)
            else: self.fromfile(filename, verbose=verbose, preview_only=preview_only)
    def _get_AtBt(self, key):
        assert(key in ['At','Bt'])
        tmp = _O.RedundantInfo.__getattribute__(self, key+'sparse')
        matrix = np.zeros((self.nAntenna + self.nUBL, len(self.crossindex)))
        for i in tmp: matrix[i[0],i[1]] = i[2]
        return sps.csr_matrix(matrix)
    def _set_AtBt(self, key, val):
        assert(key in ['At','Bt'])
        nonzeros = np.array(val.nonzero()).transpose()
        self.__setattr__(key+'sparse', np.array([[i,j,val[i,j]] for i,j in nonzeros], dtype=np.int32))
    def __getattribute__(self, key):
        if key in ['At','Bt']: return self._get_AtBt(key)
        else: return _O.RedundantInfo.__getattribute__(self, key)
    def __setattr__(self, key, val):
        if key in ['At','Bt']: return self._set_AtBt(key, val)
        else: return _O.RedundantInfo.__setattr__(self, key, val)
    #def keys(self): return [k for k in dir(self) if not k.startswith('_')] # XXX should exclude functions
    def __getitem__(self,k): return self.__getattribute__(k)
    def __setitem__(self,k,val): return self.__setattr__(k,val)
    def fromfile(self, filename, verbose=False, preview_only=False): # XXX what is preview?
        '''Initialize from (binary) file.'''
        if verbose: print 'Reading redundant info from %s' % filename
        datachunk = np.fromfile(filename)
        markerindex = np.where(datachunk == MARKER)[0]
        # XXX uneven array
        d = np.array([np.array(datachunk[markerindex[i]+1:markerindex[i+1]]) for i in range(len(markerindex)-1)])
        self.from_array(d, verbose=verbose, preview_only=preview_only) # XXX do i need to store preview return case?
        self.update()
        if verbose: print "done. nAntenna,nUBL,nBaseline = %i,%i,%i" % (len(self.subsetant),self.nUBL,self.nBaseline)
    def fromfile_txt(self, filename, verbose=False):
        '''Initialize from (txt) file.'''
        if verbose: print 'Reading redundant info from %s' % filename
        d = np.array([np.array(map(float, line.split())) for line in open(filename)])
        #assert(len(d) >= len(self.keys)) # did we get the expected number of rows?  XXX is this check necessary?
        self.from_array(d, verbose=verbose)
        self.update(txtmode=True)
        if verbose: print "done. nAntenna,nUBL,nBaseline = %i,%i,%i" % (len(self.subsetant),self.nUBL,self.nBaseline)
    def from_array(self, d, verbose=False, preview_only=False): 
        '''Initialize fields from data contained in 2D float array used to store data to file.'''
        # XXX from_array and to_array do not match, need to change that, but this affects fromfile & fromfile_txt
        self.nAntenna = int(d[0][0]) # XXX did we mean d[0,0]?
        self.nUBL = int(d[1][0]) # XXX
        self.nBaseline = int(d[2][0]) # XXX
        self.subsetant = d[3].astype(np.int32) # index of good antennas
        self.antloc = d[4].reshape((self.nAntenna,3)).astype(np.float32)
        self.subsetbl = d[5].astype(np.int32) # index of good bls (+autos) within all bls
        self.ubl = d[6].reshape((self.nUBL,3)).astype(np.float32) # unique bl vectors
        self.bltoubl = d[7].astype(np.int32) # cross bl number to ubl index
        self.reversed = d[8].astype(np.int32) # cross only bl if reversed -1, else 1
        self.reversedauto = d[9].astype(np.int32) # XXX check comment: index of good autos within all bls
        self.autoindex = d[10].astype(np.int32) # index of auto bls among good bls
        self.crossindex = d[11].astype(np.int32) # index of cross bls among good bls
        ncross = len(self.crossindex)
        # XXX maybe add this as a function
        if preview_only: return ncross - self.nUBL - self.nAntenna + 2 # XXX return value here, normally not returning anything
        self.bl2d = d[12].reshape(self.nBaseline,2).astype(np.int32) # 1d bl index to (i,j) antenna pair
        self.ublcount = d[13].astype(np.int32) # for each ubl, number of corresponding good cross bls
        self.ublindex = d[14].reshape(ncross,3).astype(np.int32) # for each ubl, the vector<int> contains (i,j,ant1,ant2,crossbl)
        self.bl1dmatrix = d[15].reshape((self.nAntenna,self.nAntenna)).astype(np.int32) #a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
        self.degenM = d[16].reshape((self.nAntenna+self.nUBL,self.nAntenna)).astype(np.float32)
        if self.nAntenna > self.threshold:
            #sparse_entries = d[16].reshape((len(d[16])/3,3))
            sparse_entries = d[17].reshape((-1,3))
            row,column,value = sparse_entries[:,0],sparse_entries[:,1],sparse_entries[:,2]
            self.At = sps.csr_matrix((value,(row,column)),shape=(ncross, self.nAntenna + self.nUBL)).T
            #sparse_entries = d[17].reshape((len(d[17])/3,3))
            sparse_entries = d[18].reshape((-1,3))
            row,column,value = sparse_entries[:,0],sparse_entries[:,1],sparse_entries[:,2]
            self.Bt = sps.csr_matrix((value,(row,column)),shape=(ncross, self.nAntenna + self.nUBL)).T
            self.AtAi = d[19].reshape((self.nAntenna + self.nUBL,self.nAntenna + self.nUBL)).astype(np.float32)
            self.BtBi = d[20].reshape((self.nAntenna + self.nUBL,self.nAntenna + self.nUBL)).astype(np.float32)
            self.totalVisibilityId = d[21].reshape(-1,2).astype(np.int32)
        else:
            # XXX why astype(int) here, but not above?
            self.At = sps.csr_matrix(d[17].reshape((ncross,self.nAntenna+self.nUBL)).astype(np.int32)).T # A matrix for logcal amplitude
            self.Bt = sps.csr_matrix(d[18].reshape((ncross,self.nAntenna+self.nUBL)).astype(np.int32)).T # B matrix for logcal phase
            self.totalVisibilityId = d[19].reshape(-1,2).astype(np.int32)
    def update(self, txtmode=False):
        '''Initialize other arrays from fundamental arrays'''
        #The sparse matrices are treated a little differently because they are not rectangular
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=DeprecationWarning)
            #if self.nAntenna <= self.threshold:
            if self.AtAi.size == 0:
                self.AtAi = la.pinv(self.At.dot(self.At.T).todense(),rcond=1e-6).astype(np.float32)#(AtA)^-1
                self.BtBi = la.pinv(self.Bt.dot(self.Bt.T).todense(),rcond=1e-6).astype(np.float32)#(BtB)^-1
    def tofile(self, filename, overwrite=False, verbose=False):
        '''XXX DOCSTRING'''
        assert(not os.path.exists(filename) or overwrite)
        if verbose: print 'Writing info to', filename
        d = self.to_array(verbose=verbose)
        f = open(filename,'wb')
        d.tofile(f)
        f.close()
        if verbose: print "Info file successfully written to", filename
    def to_array(self, verbose=False):
        # XXX from_array and to_array do not match, need to change that, but this affects fromfile & fromfile_txt
        d = ['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','At','Bt','AtAi','BtBi','totalVisibilityId']
        if self.nAntenna <= self.threshold: d = d[:-3] + d[-1:]
        def fmt(k):
            if k in ['At','Bt']: 
                sk = self[k].T
                if self.nAntenna > self.threshold:
                    row,col = sk.nonzero()
                    return np.vstack([row,col,sk[row,col]]).T
                else: return sk.todense()
            else: return self[k]
        d = [fmt(k) for k in d]
        d = [[MARKER]]+[k for i in zip(d,[[MARKER]]*len(d)) for k in i]
        return np.concatenate([np.asarray(k).flatten() for k in d])
    def compare(self, info, verbose=False, tol=1e-5):
        '''compare with another RedundantInfo, output True if they are the same and False if they are different'''
        try:
            floatkeys = float_infokeys#['antloc','ubl','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
            intkeys = ['nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix']
            infomatrices=['At','Bt']
            specialkeys = ['ublindex']
            allkeys = floatkeys + intkeys + infomatrices + specialkeys#['antloc','ubl','nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB','A','B','At','Bt']
            diff = []
            for key in floatkeys:
                if verbose: print key, 'checking'
                #ok = np.allclose(self[key], info[key], tol)
                ok = round(la.norm(np.array(self[key])-np.array(info[key]))/tol) == 0
                if verbose: print key, ok
                if not ok: return False
                #try: diff.append(round(la.norm(np.array(self[key])-np.array(info[key]))/tol)==0)
                #except: diff.append(False)
            for key in intkeys:
                if verbose: print key, 'checking'
                ok = la.norm(np.array(self[key])-np.array(info[key])) == 0
                if verbose: print key, ok
                if not ok: return False
                #try: diff.append(la.norm(np.array(self[key])-np.array(info[key]))==0)
                #except: diff.append(False)
            for key in infomatrices:
                if verbose: print key, 'checking'
                # XXX made a switch here. ok?
                #ok = la.norm(np.array(self[key])-np.array(info[key])) == 0
                #ok = la.norm((self[key]-info[key]).todense()) == 0
                ok = np.all(self[key].todense() == info[key].todense())
                if verbose: print key, ok
                if not ok: return False
                #try: diff.append(la.norm((self[key]-info[key]).todense())==0)
                #except: diff.append(False)
            diff.append(True)
            
            # XXX changed ublindex representation.  how much depends on this?
            #try:
            #    for i,j in zip(self['ublindex'],info['ublindex']):
            #        diff[-1] = diff[-1] and (la.norm(np.array(i) - np.array(j))==0)
            #except: diff[-1] = False
            bool = True
            for i in diff: bool &= i
            #print the first key found different (this will only trigger when the two info's have the same shape, so probably not very useful)
            if verbose and not bool:
                for i in range(len(diff)):
                    if not diff[i]: print allkeys[i]
            return bool
        except(ValueError):
            print "info doesn't have the same shape"
            return False
    def get_xy_AB(self):
        '''return xyA, xyB, yxA, yxB for logcal cross polarizations'''
        na = self.nAntenna
        nu = self.nUBL
        A = self.At.T.todense()
        B = self.Bt.T.todense()
        bl2dcross = self.bl2d[self.crossindex]
        #print na, nu, B.shape wesdcxaz
        xyA = np.zeros((len(self.crossindex), 2*na+nu), dtype='int8')
        yxA = np.zeros_like(xyA)
        xyB = np.zeros_like(xyA)
        yxB = np.zeros_like(xyA)
        xyA[:, 2*na:] = A[:, na:]
        xyB[:, 2*na:] = B[:, na:]
        for i in range(len(xyA)):
            xyA[i, bl2dcross[i,0]] = A[i, bl2dcross[i,0]]
            xyA[i, na + bl2dcross[i,1]] = A[i, bl2dcross[i,1]]
            xyB[i, bl2dcross[i,0]] = B[i, bl2dcross[i,0]]
            xyB[i, na + bl2dcross[i,1]] = B[i, bl2dcross[i,1]]
        yxA[:, :na] = xyA[:, na:2*na]
        yxA[:, na:2*na] = xyA[:, :na]
        yxA[:, 2*na:] = xyA[:, 2*na:]
        yxB[:, :na] = xyB[:, na:2*na]
        yxB[:, na:2*na] = xyB[:, :na]
        yxB[:, 2*na:] = xyB[:, 2*na:]
        return xyA, xyB, yxA, yxB
    
