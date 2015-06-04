import _omnical as _O
import numpy as np, numpy.linalg as la
import warnings, os, time
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy.sparse as sps

KEYS = [
    'nAntenna', # number of usable ants (not total number)
    'nUBL', # number of unique bls, matches first dim of ubl/ublcount/ublindex
    'nBaseline', # number of bls, matches first dim of bltoubl/bl2d, now python only
    'subsetant', # (nAntenna,) antenna numbers used; index i corresponds to antenna number ai, now python only
    'antloc', # (nAntenna,3) float,  idealized antpos from which redundancy is determined XXX not sure of lin/log cal need this.  if not, could restrict this to ArrayInfo and remove from RedundantInfo
    'subsetbl', # (nBaseline,) for each bl in bl2d, the index in totalVisibilityId; now python only
    'ubl', # (nUBL,3) float, sep vector for each unique baseline i think not necessary for lin/log cal, now python only
    'bltoubl', # (nBaseline,) for each bl in bl2d, the index of corresponding unique bl in ubl/ublcount/ublindex
    'reversed', # for each entry in crossindex, 1 if baseline is flipped wrt corresponding ubl, otherwise -1
    'reversedauto', # XXX to read old files
    'autoindex', # XXX to read old files
    'crossindex', # indices in bl2d of crosses XXX if we mandate no autos, then crossindex not necessary
    'bl2d', # (nBaseline,2) the i,j indices of ants in subsetant for each bl
    'ublcount', # (nUBL,) number of bls contributing to each ubl XXX can determine this from ublindex
    'ublindex', # (nUBL, ublcount[i], 3) ant1,ant2,blindex for each bl contributing to each ubl
    'bl1dmatrix', # (nAntenna,nAntenna) for each i,j antenna pair, the index of where that bl appears in crossindex
    'degenM', # (nAntenna+nUBL,nAntenna)
    'At', # (ncross,nAntenna+nUBL), sparse
    'Bt', # (ncross,nAntenna+nUBL), sparse
    'AtAi', # precomputed matrix
    'BtBi', # precomputed matrix
    'totalVisibilityId', # (all_baselines, 2) ai,aj antenna numbers for every possible bl; defines order of data to be loaded into omnical solver XXX if totalVisibilityId only holds good data, then subsetbl becomes pointless
]

# XXX idea for better interface: provide list of bls grouped by ubl & whether they need to be flipped
# this list must already omit bad bls and ants
# determine subsetant from this list & internally take care of mapping indices to antennas.
# then do major cleanup of c code

# XXX all this meta stuff about "info" almost assuredly means info needs to be a class
#infokeys = ['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B','At','Bt','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
infokeys = ['nAntenna','nUBL','antloc','bltoubl','reversed','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B','At','Bt','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
infokeys_optional = ['totalVisibilityId']
#binaryinfokeys=['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','At','Bt']
binaryinfokeys=['nAntenna','nUBL','antloc','bltoubl','reversed','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','At','Bt']
cal_name = {0: "Lincal", 1: "Logcal"}

int_infokeys = ['nAntenna','nUBL']
#intarray_infokeys = ['subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','A','B','At','Bt']
intarray_infokeys = ['bltoubl','reversed','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','A','B','At','Bt']
intarray_infokeys_optional = ['totalVisibilityId']

#float_infokeys = ['antloc','ubl','degenM','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
float_infokeys = ['degenM','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']

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
        self.totalVisibilityId = np.zeros_like(self.bl2d) # XXX placeholder for now
    def _get_AtBt(self, key):
        '''for convenience of multiplication in update()'''
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
    def to_npz(self, filename):
        def fmt(k):
            if k in ['At','Bt']: return _O.RedundantInfo.__getattribute__(self,k+'sparse')
            else: return self[k]
        d = {}
        for k in KEYS: d[k] = fmt(k)
        np.savez(filename, **d)
    def from_npz(self, filename):
        npz = np.load(filename)
        def fmt(npz,k):
            if k in ['nAntenna','nUBL','nBaseline']: self[k] = int(npz[k])
            elif k in ['At','Bt']: _O.RedundantInfo.__setattr__(self,k+'sparse', npz[k])
            else: self[k] = npz[k]
        for k in KEYS: fmt(npz,k)
        self.update()
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
        # XXX maybe not all of these variables should be exposed (i.e. some could be C only)
        # XXX at the least, should validate array dimensions in C wrapper
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
        #self.crossindex = d[11].astype(np.int32) # index of cross bls among good bls
        crossindex = d[11].astype(np.int32) # index of cross bls among good bls
        self.crossindex = np.arange(len(crossindex), dtype=np.int32)
        ncross = len(self.crossindex)
        # XXX maybe add this as a function
        if preview_only: return ncross - self.nUBL - self.nAntenna + 2 # XXX return value here, normally not returning anything
        #self.bl2d = d[12].reshape(self.nBaseline,2).astype(np.int32) # 1d bl index to (i,j) antenna pair
        bl2d = d[12].reshape(self.nBaseline,2).astype(np.int32) # 1d bl index to (i,j) antenna pair
        self.bl2d = bl2d[crossindex]
        self.ublcount = d[13].astype(np.int32) # for each ubl, number of corresponding good cross bls
        #self.ublindex = d[14].reshape(ncross,3).astype(np.int32) # for each ubl, the vector<int> contains (i,j,ant1,ant2,crossbl)
        ublindex = d[14].reshape(ncross,3).astype(np.int32) # for each ubl, the vector<int> contains (i,j,ant1,ant2,crossbl)
        newind = np.arange(self.nBaseline)[crossindex] = self.crossindex
        ublindex[2] = newind[ublindex[2]]
        self.ublindex = ublindex
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
        d = KEYS
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
    def init_from_redundancies(self, reds, antpos):
        '''Initialize RedundantInfo from a list where each entry is a group of redundant baselines.
        each baseline is a (i,j) tuple, where i,j are antenna indices.  To ensure baselines are
        oriented to be redundant, it may be necessary to have i > j.  If this is the case, then
        when calibrating visibilities listed as j,i data will have to be conjugated.'''
        ants = {}
        for ubl_gp in reds:
            for (i,j) in ubl_gp: ants[i] = ants[j] = None
        self.subsetant = np.array(ants.keys(), dtype=np.int32)
        ant2ind = {}
        for i,ant in enumerate(self.subsetant): ant2ind[ant] = i
        self.nAntenna = self.subsetant.size
        self.nUBL = len(reds)
        #bl2d = np.array([(ant2ind[i],ant2ind[j],u,i<j) for u,ubl_gp in enumerate(reds) for i,j in ubl_gp], dtype=np.int32)
        bl2d = np.array([(ant2ind[min(i,j)],ant2ind[max(i,j)],u,i<j) for u,ubl_gp in enumerate(reds) for i,j in ubl_gp], dtype=np.int32)
        self.bl2d = bl2d[:,:2]
        self.nBaseline = bl2d.shape[0]
        self.bltoubl = bl2d[:,2]
        #self.reversed = np.where(bl2d[:,3] == 0, 1, -1).astype(np.int32)
        self.reversed = np.where(bl2d[:,3] == 0, -1, 1).astype(np.int32)
        self.crossindex = np.arange(self.nBaseline, dtype=np.int32) # XXX mandating no autos here
        # XXX subsetbl
        self.reversedauto = np.ones_like(self.reversed)
        self.autoindex = np.ones_like(self.crossindex)
        self.subsetbl = np.arange(self.nBaseline, dtype=np.int32) # XXX mandating visibilities provided in same order
        # remvoe above eventually
        self.ublcount = np.array([len(ubl_gp) for ubl_gp in reds], dtype=np.int32)
        bl2d[:,2] = np.arange(self.nBaseline)
        self.ublindex = bl2d[:,:3]
        bl1dmatrix = (2**31-1) * np.ones((self.nAntenna,self.nAntenna),dtype=np.int32)
        for n,(i,j) in enumerate(self.bl2d): bl1dmatrix[i,j], bl1dmatrix[j,i] = n,n
        self.bl1dmatrix = bl1dmatrix
        #A: A matrix for logcal amplitude
        A,B = np.zeros((self.nBaseline,self.nAntenna+self.nUBL)), np.zeros((self.nBaseline,self.nAntenna+self.nUBL))
        for n,(i,j) in enumerate(self.bl2d):
            A[n,i], A[n,j], A[n,self.nAntenna+self.bltoubl[n]] = 1,1,1
            B[n,i], B[n,j], B[n,self.nAntenna+self.bltoubl[n]] = -self.reversed[n],self.reversed[n],1
        self.At, self.Bt = sps.csr_matrix(A).T, sps.csr_matrix(B).T
        # XXX nothing up to this point requires antloc
        self.antloc = antpos.take(self.subsetant, axis=0).astype(np.float32) # XXX check this
        self.ubl = np.array([np.mean([antpos[j]-antpos[i] for i,j in ubl_gp],axis=0) for ubl_gp in reds], dtype=np.float32) # XXX check this
        # XXX would like to understand better what is happening here
        a = np.array([np.append(ai,1) for ai in self.antloc], dtype=np.float32)
        d = np.array([np.append(ubli,0) for ubli in self.ubl], dtype=np.float32)
        m1 = -a.dot(la.pinv(a.T.dot(a))).dot(a.T)
        m2 = d.dot(la.pinv(a.T.dot(a))).dot(a.T)
        self.degenM = np.append(m1,m2,axis=0)
        self.update()
    def list_redundancies(self):
        '''After initialization, return redundancies in the same format used in init_from_redundancies.'''
        # XXX broken
        #bls = self.bl2d[self.crossindex,:2]
        #bls = [(i,j) if self.reversed[cnt] < 0 else (j,i) for cnt,(i,j) in enumerate(bls)]
        reds = []
        x = 0
        for y in self.ublcount:
            #reds.append([(self.subsetant[i],self.subsetant[j]) for i,j in self.ublindex[x:x+y,:2]])
            reds.append([(self.subsetant[i],self.subsetant[j]) if self.reversed[k] == 1 else (self.subsetant[j],self.subsetant[i]) for i,j,k in self.ublindex[x:x+y]])
            x += y
        return reds
    #def from_arrayinfo(self, arrayinfoPath=None, verbose=False, badAntenna=[], badUBLpair=[], tol=1e-6):
    #    '''Use provided antenna locations (in arrayinfoPath) to derive redundancy equations'''
    #    if arrayinfoPath is not None: self.read_arrayinfo(arrayinfoPath)
    #    # exclude bad antennas
    #    self['subsetant'] = subsetant = np.array([i for i in xrange(antennaLocation.shape[0]) 
    #            if i not in badAntenna], dtype=np.int32)
    #    self['nAntenna'] = nAntenna = len(subsetant) # XXX maybe have C api automatically infer this
    #    self['antloc'] = antloc = np.array([antennaLocation[i] for i in subsetant], dtype=np.float32)
    #    #delete the bad ubl's
    #    badUBL = {}
    #    def dis(a1,a2): return np.linalg.norm(a1-a2)
    #    for a1,a2 in badUBLpair:
    #        bl = antennaLocation[a1] - antennaLocation[a2]
    #        for i,ubl in enumerate(ublall):
    #            if dis(bl,ubl) < tol or dis(bl,-ubl) < tol: badUBL[i] = None
    #    ubl2goodubl = {}
    #    def f(i,u):
    #        ubl2goodubl[i] = len(ubl2goodubl)
    #        return u
    #    self['ubl'] = ubl = np.array([f(i,u) for i,u in enumerate(ublall) if not badUBL.has_key(i)], dtype=np.float32)
    #    for k in badUBL: ubl2goodubl[k] = -1
    #    self['nUBL'] = nUBL = ubl.shape[0] # XXX maybe have C api automatically infer this
    #    badubl = [ublall[i] for i in badUBL]
    #    #find nBaseline (include auto bls) and subsetbl
    #    #bl2d:  from 1d bl index to a pair of antenna numbers
    #    bl2d = [] # XXX cleaner way to do this?
    #    for i,ai in enumerate(antloc):
    #        for j,aj in enumerate(antloc[:i+1]):
    #            blij = ai - aj
    #            flag = False
    #            for bl in badubl:
    #                if dis(blij,bl) < tol or dis(blij,-bl) < tol:
    #                    flag = True
    #                    break
    #            if not flag: bl2d.append((i,j))
    #    # exclude pairs that are not in totalVisibilityId
    #    tmp = []
    #    for p in bl2d:
    #        bl = (subsetant[p[0]],subsetant[p[1]])
    #        if totalVisibilityId_dic.has_key(bl): tmp.append(p)
    #        elif totalVisibilityId_dic.has_key(bl[::-1]): tmp.append(p[::-1])
    #    self['bl2d'] = bl2d = np.array(tmp, dtype=np.int32)
    #    self['nBaseline'] = len(bl2d) # XXX maybe have C api infer this
    #    # from a pair of good antenna index to bl index
    #    self['subsetbl'] = np.array([self.get_baseline([subsetant[bl[0]],subsetant[bl[1]]]) 
    #            for bl in bl2d], dtype=np.int32)
    #    #bltoubl: cross bl number to ubl index
    #    def findublindex(p1,p2):
    #        a1,a2 = subsetant[p1],subsetant[p2]
    #        if (a1,a2) in self.totalVisibilityUBL: return ubl2goodubl[self.totalVisibilityUBL[(a1,a2)]]
    #    self['bltoubl'] = bltoubl = np.array([findublindex(*p) for p in bl2d if p[0] != p[1]], dtype=np.int32)
    #    #reversed:   cross only bl if reversed -1, otherwise 1
    #    crosspair = [p for p in bl2d if p[0] != p[1]]
    #    reverse = []
    #    for k,cpk in enumerate(crosspair):
    #        bl = antloc[cpk[0]] - antloc[cpk[1]]
    #        if dis(bl,ubl[bltoubl[k]]) < tol: reverse.append(-1)
    #        elif dis(bl,-ubl[bltoubl[k]]) < tol: reverse.append(1)
    #        else : raise ValueError('bltoubl[%d] points to wrong ubl index' % (k))
    #    self['reversed'] = np.array(reverse, dtype=np.int32)
    #    # autoindex: index of auto bls among good bls
    #    self['autoindex'] = autoindex = np.array([i for i,p in enumerate(bl2d) if p[0] == p[1]], dtype=np.int32)
    #    # crossindex: index of cross bls among good bls
    #    self['crossindex'] = crossindex = np.array([i for i,p in enumerate(bl2d) if p[0] != p[1]], dtype=np.int32)
    #    # reversedauto: the index of good bls (auto included) in all bls
    #    reversedauto = np.arange(self['nBaseline'], dtype=np.int32)
    #    for i in autoindex: reversedauto[i] = 1
    #    for i,c in enumerate(crossindex): reversedauto[c] = reverse[i]
    #    self['reversedauto'] = reversedauto
    #    #ublcount:  for each ubl, the number of good cross bls corresponding to it
    #    cnt = {}
    #    for bl in bltoubl: cnt[bl] = cnt.get(bl,0) + 1
    #    self['ublcount'] = np.array([cnt[i] for i in range(nUBL)], dtype=np.int32)
    #    #ublindex:  //for each ubl, the vector<int> contains (ant1, ant2, crossbl)
    #    cnt = {}
    #    for i,(a1,a2) in enumerate(crosspair): cnt[bltoubl[i]] = cnt.get(bltoubl[i],[]) + [[a1,a2,i]]
    #    self['ublindex'] = ublindex = np.concatenate([np.array(cnt[i],dtype=np.int32) for i in range(nUBL)])
    #    #bl1dmatrix: a symmetric matrix where col/row numbers index ants and entries are bl index (no auto corr)
    #    # XXX don't like 2**31-1.  whence this number?
    #    bl1dmatrix = (2**31-1) * np.ones((nAntenna,nAntenna),dtype=np.int32)
    #    for i,cp in enumerate(crosspair): bl1dmatrix[cp[1],cp[0]], bl1dmatrix[cp[0],cp[1]] = i,i
    #    self['bl1dmatrix'] = bl1dmatrix
    #    #degenM:
    #    a = np.array([np.append(ai,1) for ai in antloc], dtype=np.float32)
    #    d = np.array([np.append(ubli,0) for ubli in ubl], dtype=np.float32)
    #    m1 = -a.dot(la.pinv(a.T.dot(a))).dot(a.T)
    #    m2 = d.dot(la.pinv(a.T.dot(a))).dot(a.T)
    #    self['degenM'] = np.append(m1,m2,axis=0)
    #    #A: A matrix for logcal amplitude
    #    A = np.zeros((len(crosspair),nAntenna+nUBL))
    #    for i,cp in enumerate(crosspair): A[i,cp[0]], A[i,cp[1]], A[i,nAntenna+bltoubl[i]] = 1,1,1
    #    self['At'] = sps.csr_matrix(A).T
    #    #B: B matrix for logcal phase
    #    B = np.zeros((len(crosspair),nAntenna+nUBL))
    #    for i,cp in enumerate(crosspair): B[i,cp[0]], B[i,cp[1]], B[i,nAntenna+bltoubl[i]] = -reverse[i],reverse[i],1
    #    self['Bt'] = sps.csr_matrix(B).T
    #    self.update()
    def compare(self, info, verbose=False, tol=1e-5):
        '''compare with another RedundantInfo, output True if they are the same and False if they are different'''
        try:
            floatkeys = float_infokeys#['antloc','ubl','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
            intkeys = ['nAntenna','nUBL','bltoubl','reversed','crossindex','bl2d','ublcount','bl1dmatrix']
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
    
