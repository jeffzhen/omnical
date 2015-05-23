import _omnical as _O
import numpy as np, numpy.linalg as la
from array import array # XXX what does array do that np.array does not?
import warnings, os
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy.sparse as sps

# XXX all this meta stuff about "info" almost assuredly means info needs to be a class
infokeys = ['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B','At','Bt','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
infokeys_optional = ['totalVisibilityId']
binaryinfokeys=['nAntenna','nUBL','nBaseline','subsetant','antloc','subsetbl','ubl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','degenM','A','B']
cal_name = {0: "Lincal", 1: "Logcal"}

int_infokeys = ['nAntenna','nUBL','nBaseline']
intarray_infokeys = ['subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','ublindex','bl1dmatrix','A','B','At','Bt']
intarray_infokeys_optional = ['totalVisibilityId']

float_infokeys = ['antloc','ubl','degenM','AtAi','BtBi']#,'AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']

def dis(a1,a2): # XXX this never used?
    '''calculate the norm of the difference of two vectors'''
    return np.linalg.norm(np.array(a1)-np.array(a2))

#  ___        _              _          _   ___       __     
# | _ \___ __| |_  _ _ _  __| |__ _ _ _| |_|_ _|_ _  / _|___ 
# |   / -_) _` | || | ' \/ _` / _` | ' \  _|| || ' \|  _/ _ \
# |_|_\___\__,_|\_,_|_||_\__,_\__,_|_||_\__|___|_||_|_| \___/

class RedundantInfo(_O.RedundantInfo):
    '''Container for metadata used by redundant calibrator, which is eventually passed into C++ routines.'''
    def __init__(self, filename=None, verbose=False, preview_only=False, threshold=128):
        _O.RedundantInfo.__init__(self)
        self.threshold = threshold # XXX move to versions
        if filename: self.fromfile(filename, verbose=verbose, preview_only=preview_only)
        # XXX don't know if we should be in the business of automatic casting.
        #elif key in int_infokeys:
        #    self.__setattr__(key, int(info[key]))
        #elif key in intarray_infokeys and key != 'ublindex':
        #    self.__setattr__(key, np.array(info[key]).astype('int32'))
        #elif key in intarray_infokeys_optional:
        #    try:
        #        self.__setattr__(key, np.array(info[key]).astype('int32'))
        #    except KeyError:
        #        pass
        #elif key in float_infokeys:
        #    self.__setattr__(key, np.array(info[key]).astype('float32'))
    def _get_AB(self, key): # XXX access to A,B is deprecated.  should remove
        assert(key in ['A','B'])
        return sps.csr_matrix(_O.RedundantInfo.__getattribute__(self, key))
    def _get_AtBt(self, key):
        assert(key in ['At','Bt'])
        tmp = _O.RedundantInfo.__getattribute__(self, key+'sparse')
        matrix = np.zeros((self.nAntenna + self.nUBL, len(self.crossindex)))
        for i in tmp: matrix[i[0],i[1]] = i[2]
        return sps.csr_matrix(matrix)
    def _set_AtBt(self, key, val):
        assert(key in ['At','Bt'])
        tmp = []
        nonzeros = np.array(val.nonzero()).transpose()
        for i,j in nonzeros:tmp += [[i, j, val[i,j]]]
        #for i in range(info[key].shape[0]):
            #for j in range(info[key].shape[1]):
                #if info[key][i,j] != 0:
                    #tmp += [[i, j, info[key][i,j]]]
        self.__setattr__(key+'sparse', np.array(tmp, dtype = 'int32'))
    def _get_ublindex(self):
        ublindex = []
        for i in self.ublindex: # XXX careful about hooking this method to __getattribute__ !
            while len(ublindex) < i[0] + 1: ublindex.append(np.zeros((1,3)))
            while len(ublindex[i[0]]) < i[1] + 1: ublindex[i[0]] = np.array(ublindex[i[0]].tolist() + [[0,0,0]])
            ublindex[i[0]][i[1]] = np.array(i[2:])
        return ublindex
    def __getattribute__(self, key):
        if key in ['A','B']: return self._get_AB(key) # XXX depricated
        elif key in ['At','Bt']: return self._get_AtBt(key)
        #elif key in ['ublindex']: # XXX
            # self._get_ublindex(self)
        else: return _O.RedundantInfo.__getattribute__(self, key)
    def __setattr__(self, key, val):
        #if key in ['A','B']: return self._set_AB(key,val) # XXX depricated
        if key in ['At','Bt']: return self._set_AtBt(key, val)
        #elif key in ['ublindex']: # XXX
            # self._get_ublindex(self)
        else: return _O.RedundantInfo.__setattr__(self, key, val)
    #def keys(self): return [k for k in dir(self) if not k.startswith('_')] # XXX should exclude functions
    def __getitem__(self,k): return self.__getattribute__(k)
    def __setitem__(self,k,val): return self.__setattr__(k,val)
    def fromfile(self, filename, verbose=False, preview_only=False): # XXX what is preview?
        '''Initialize from (binary) file.'''
        if verbose:
            print 'Reading redundant info from %s' % filename
            time = time.time()
        f = open(filename)
        farray = array('d')
        farray.fromstring(f.read())
        datachunk = np.array(farray)
        marker = 9999999
        markerindex = np.where(datachunk == marker)[0]
        # XXX uneven array
        d = np.array([np.array(datachunk[markerindex[i]+1:markerindex[i+1]]) for i in range(len(markerindex)-1)])
        self.from_array(d, verbose=verbose, preview_only=preview_only) # XXX do i need to store preview return case?
        self.update()
        if verbose:
            print "done. nAntenna, nUBL, nBaseline = %i, %i, %i. Time taken: %f minutes." % (len(self.subsetant),self.nUBL,self.nBaseline,(time.time()-timer)/60.)
    def fromfile_txt(self, filename, verbose=False):
        '''Initialize from (txt) file.'''
        if verbose: 
            print 'Reading redundant info from %s' % filename
            time = time.time()
        f = open(filename)
        d = np.array([np.array(map(float, line.split())) for line in f])
        #assert(len(d) >= len(self.keys)) # did we get the expected number of rows?  XXX is this check necessary?
        self.from_array(d, verbose=verbose)
        self.update(txtmode=True)
        if verbose:
            print "done. nAntenna, nUBL, nBaseline = %i, %i, %i. Time taken: %f minutes." % (len(self.subsetant),self.nUBL,self.nBaseline,(time.time()-timer)/60.)
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
        # XXX i think this could be done more elegantly
        try: # XXX this is hacky for backward compatibility
            tmp = d[14].reshape(ncross,3).astype(np.int32)
            cnt = 0
            ublindex = [] # for each ubl, the vector<int> contains (i,j,ant1,ant2,crossbl)
            for i in range(self.nUBL):
                for j in range(self.ublcount[i]):
                    ublindex.append([i,j,tmp[cnt][0],tmp[cnt][1],tmp[cnt][2]])
        except(ValueError): # means this is the new ublindex that already has i,j in it
            ublindex = d[14].reshape(ncross,5).astype(np.int32)
        self.ublindex = np.asarray(ublindex, dtype=np.int32)
        self.bl1dmatrix = d[15].reshape((self.nAntenna,self.nAntenna)).astype(np.int32) #a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
        self.degenM = d[16].reshape((self.nAntenna+self.nUBL,self.nAntenna)).astype(np.float32)
        if self.nAntenna > self.threshold:
            #sparse_entries = d[16].reshape((len(d[16])/3,3))
            sparse_entries = d[17].reshape((-1,3))
            row,column,value = sparse_entries[:,0],sparse_entries[:,1],sparse_entries[:,2]
            # XXX A and B are never accessed, only At and Bt
            self.A = sps.csr_matrix((value,(row,column)),shape=(ncross, self.nAntenna + self.nUBL))
            #sparse_entries = d[17].reshape((len(d[17])/3,3))
            sparse_entries = d[18].reshape((-1,3))
            row,column,value = sparse_entries[:,0],sparse_entries[:,1],sparse_entries[:,2]
            self.B = sps.csr_matrix((value,(row,column)),shape=(ncross, self.nAntenna + self.nUBL))
            self.AtAi = d[19].reshape((self.nAntenna + self.nUBL,self.nAntenna + self.nUBL))
            self.BtBi = d[20].reshape((self.nAntenna + self.nUBL,self.nAntenna + self.nUBL))
            self.totalVisibilityId = d[21].reshape(-1,2).astype(np.int32)
        else:
            # XXX why astype(int) here, but not above?
            self.A = sps.csr_matrix(d[17].reshape((ncross,self.nAntenna+self.nUBL)).astype(np.int32)) # A matrix for logcal amplitude
            self.B = sps.csr_matrix(d[18].reshape((ncross,self.nAntenna+self.nUBL)).astype(np.int32)) # B matrix for logcal phase
            self.totalVisibilityId = d[19].reshape(-1,2).astype(np.int32)
    def update(self, txtmode=False):
        '''Initialize other arrays from fundamental arrays'''
        #The sparse matrices are treated a little differently because they are not rectangular
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",category=DeprecationWarning)
            self.At = self.A.T # do we really need this? yes, b/c A is never used
            self.Bt = self.B.T 
            if self.nAntenna <= self.threshold:
                self.AtAi = la.pinv(self.At.dot(self.A).todense(),rcond=1e-6).astype(np.float32)#(AtA)^-1
                self.BtBi = la.pinv(self.Bt.dot(self.B).todense(),rcond=1e-6).astype(np.float32)#(BtB)^-1
            if txtmode: # XXX from here on, only in txt, probably not useful?
                self.AtAiAt = self.AtAi.dot(self.At.todense())#(AtA)^-1At
                self.BtBiBt = self.BtBi.dot(self.Bt.todense())#(BtB)^-1Bt
                self.PA = self.A.dot(self.AtAiAt)#A(AtA)^-1At
                self.PB = self.B.dot(self.BtBiBt)#B(BtB)^-1Bt
                ncross = len(self.crossindex)
                self.ImPA = sps.identity(ncross) - self.PA#I-PA
                self.ImPB = sps.identity(ncross) - self.PB#I-PB
    def tofile(self, filename, overwrite=False, verbose=False):
        '''XXX DOCSTRING'''
        assert(not os.path.exists(filename) or overwrite)
        if verbose:
            timer = time.time()
            print 'Writing info to', filename
        d = self.to_array(verbose=verbose)
        f = open(filename,'wb')
        d.tofile(f)
        f.close()
        if verbose:
            print "Info file successfully written to %s. Time taken: %f minutes."%(filename, (time.time()-timer)/60.)
    def to_array(self, verbose=False):
        # XXX from_array and to_array do not match, need to change that, but this affects fromfile & fromfile_txt
        binaryinfokeysnew = binaryinfokeys[:]
        if self.nAntenna > self.threshold: binaryinfokeysnew.extend(['AtAi','BtBi'])
        binaryinfokeysnew.extend(['totalVisibilityId']) # XXX propose that when writing, can break backward compatibility
        #if 'totalVisibilityId' in info.keys(): binaryinfokeysnew.extend(['totalVisibilityId'])
        #else: print "warning: info doesn't have the key 'totalVisibilityId'"
        marker = 9999999 # XXX what is this magic number for?
        #XXX need to do the rest more elegantly
        datachunk = [0 for i in range(len(binaryinfokeysnew)+1)]
        count = 0
        datachunk[count] = np.array([marker])         #start with a marker
        count += 1
        if verbose: print "appending",
        for key in binaryinfokeysnew:
            if key in ['antloc', 'ubl','degenM', 'AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB','totalVisibilityId']:  #'antloc',
                add = np.append(np.array(self[key]).flatten(),[marker])
                datachunk[count] = add
                count += 1
                if verbose: print key,
            elif key == 'ublindex':
                add = np.append(np.vstack(self[key]).flatten(),[marker])
                datachunk[count] = add
                count += 1
                if verbose: print key,
            elif key in ['A','B']: # XXX do we need to write A,B to file?
                if self.nAntenna > self.threshold:
                    row = self[key].nonzero()[0]           #row index of non zero entries
                    column = self[key].nonzero()[1]        #column index of non zero entries
                    nonzero = np.transpose(np.array([row,column]))       #a list of non zero entries
                    temp = np.array([np.array([row[i],column[i],self[key][nonzero[i,0],nonzero[i,1]]]) for i in range(len(nonzero))])
                    add = np.append(np.array(temp.flatten()),[marker])
                    datachunk[count] = add
                else:
                    add = np.append(np.array(self[key].todense().flatten()).flatten(),[marker])
                    datachunk[count] = add
                count += 1
                if verbose: print key,
            else:
                add = np.append(np.array(self[key]).flatten(),[marker])
                datachunk[count] = add
                count += 1
                if verbose: print key,
        if verbose: print ""
        datachunkarray = array('d',np.concatenate(tuple(datachunk)))
        return datachunkarray
    # XXX part of Info class?
    def compare(self, info, verbose=False, tol=1e-5):
        '''compare with another RedundantInfo, output True if they are the same and False if they are different'''
        try:
            floatkeys = float_infokeys#['antloc','ubl','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB']
            intkeys = ['nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix']
            #infomatrices=['A','B','At','Bt']
            infomatrices=['At','Bt'] # XXX access to A,B depreciated, right?
            specialkeys = ['ublindex']
            allkeys = floatkeys + intkeys + infomatrices + specialkeys#['antloc','ubl','nAntenna','nUBL','nBaseline','subsetant','subsetbl','bltoubl','reversed','reversedauto','autoindex','crossindex','bl2d','ublcount','bl1dmatrix','AtAi','BtBi','AtAiAt','BtBiBt','PA','PB','ImPA','ImPB','A','B','At','Bt']
            diff = []
            for key in floatkeys:
                if verbose: print key, 'checking'
                #ok = np.allclose(self[key], info[key], tol) # XXX segfaults
                ok = round(la.norm(np.array(self[key])-np.array(info[key]))/tol) == 0 # XXX also segfaults
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
        A = self.A.todense()
        B = self.B.todense()
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
    def remove_one_antenna(self,badant):
        '''XXX DOCSTRING'''
        # XXX this function is never called anywhere, as far as I can tell
        #nAntenna and antloc
        nAntenna = self.nAntenna-1
        badindex = list(self.subsetant).index(badant) # index of badant in prev subsetant
        subsetant = list(self.subsetant)[:]
        subsetant.pop(badindex) #delete the bad antenna from subsetant
        antloc = np.delete(np.array(self.antloc),badindex,0) #delete the bad antenna from antloc
        #ubl and nUBL
        index = 0 #to keep track of the index of ubl the loop is at
        deletelist = []
        for ubl in self.ublindex:
            if len(ubl) > 1: index += 1
            elif ubl[0,0] == subsetant[badant] or ubl[0,1] == subsetant[badant] :
                deletelist.append(index)
                index += 1
        ubl = self.ubl[:]
        ubl = np.array([ubl[i] for i in range(len(ubl)) if i not in deletelist])
        nUBL=len(ubl);
        #subsetbl and nBaseline
        goodpairs_old = [i[::-1] for i in self.bl2d]    #the old goodpairs
        #the index of goodpairs that doesn't contain the bad antenna
        goodpairs_index = [i for i in range(len(goodpairs_old)) if badindex not in goodpairs_old[i]]
        #the new goodpairs with antenna number (all antenna)
        temp = np.array([goodpairs_old[i] for i in goodpairs_index])
        goodpairs = np.zeros(temp.shape)
        for i in range(len(temp)):
            for j in range(len(temp[i])):
                if temp[i,j] > badindex: goodpairs[i,j] = temp[i,j]-1
                else: goodpairs[i,j] = temp[i,j]
        subsetbl = [self.subsetbl[i] for i in goodpairs_index]  #the new subsetbl
        nBaseline = len(subsetbl)
        counter = 0
        ubl_old2new = np.zeros([len(self.ubl)], dtype='int') #from old ubl index to new ubl index
        for i in range(len(self.ubl)):
            if i in deletelist: ubl_old2new[i] = counter
            else:
                ubl_old2new[i] = counter
                counter += 1
        bltoubl = []
        for i in range(len(self.crossindex)):
            #get the pair of antenna from each crossindex
            pair = [self.subsetant[index] for index in self.bl2d[self.crossindex[i]]]
            #append the new ubl index that doesn't have the bad antenna
            if not badant in pair: bltoubl.append(ubl_old2new[self.bltoubl[i]]) 
        bltoubl = np.array(bltoubl)
        #################################################################################
        #reversed:   cross only bl if reversed -1, otherwise 1
        crosspair_old = []
        for p in goodpairs_old:
            if p[0] != p[1]: crosspair_old.append(p)
        goodcross = []
        for i in range(len(crosspair_old)):
            if badindex not in crosspair_old[i]: goodcross.append(i)
        crosspair=[]
        for p in goodpairs:
            if p[0] != p[1]: crosspair.append(p)
        reverse=[self.reversed[i] for i in goodcross]
        ######################################################################################
        #reversedauto: the index of good baselines (auto included) in all baselines
        #autoindex: index of auto bls among good bls
        #crossindex: index of cross bls among good bls
        #ncross
        reversedauto = range(len(goodpairs))
        #find the autoindex and crossindex in goodpairs
        autoindex = []
        crossindex = []
        for i in range(len(goodpairs)):
            if goodpairs[i][0] == goodpairs[i][1]: autoindex.append(i)
            else: crossindex.append(i)
        for i in autoindex: reversedauto[i] = 1
        for i in range(len(crossindex)): reversedauto[crossindex[i]] = reverse[i]
        reversedauto = np.array(reversedauto)
        autoindex = np.array(autoindex)
        crossindex = np.array(crossindex)
        ncross = len(crossindex)
        ###################################################
        #bl2d:  from 1d bl index to a pair of antenna numbers
        bl2d = []
        for pair in goodpairs: bl2d.append(pair[::-1])
        bl2d = np.array(bl2d)
        ###################################################
        #ublcount:  for each ubl, the number of good cross bls corresponding to it
        countdict = {}
        for bl in bltoubl: countdict[bl] = 0
        for bl in bltoubl: countdict[bl] += 1
        ublcount=[]
        for i in range(nUBL): ublcount.append(countdict[i])
        ublcount = np.array(ublcount)
        ####################################################################################
        #ublindex:  //for each ubl, the vector<int> contains (ant1, ant2, crossbl)
        countdict = {}
        for bl in bltoubl: countdict[bl] = []
        for i in range(len(crosspair)):
            ant1 = crosspair[i][1]
            ant2 = crosspair[i][0]
            countdict[bltoubl[i]].append([ant1,ant2,i])
        ublindex=[]
        for i in range(nUBL): ublindex.append(countdict[i])
        #turn each list in ublindex into np array
        for i in range(len(ublindex)): ublindex[i]=np.array(ublindex[i])
        ublindex = np.array(ublindex)
        ###############################################################################
        #bl1dmatrix: a symmetric matrix where col/row numbers are antenna indices and entries are 1d baseline index not counting auto corr
        #I suppose 99999 for bad and auto baselines?
        bl1dmatrix = 99999*np.ones([nAntenna,nAntenna],dtype='int16') # XXX don't like magic numbers
        for i in range(len(crosspair)):
            bl1dmatrix[crosspair[i][1]][crosspair[i][0]] = i
            bl1dmatrix[crosspair[i][0]][crosspair[i][1]] = i
        ####################################################################################3
        #degenM:
        a = []
        for i in range(len(antloc)): a.append(np.append(antloc[i],1))
        a = np.array(a)
        d = []
        for i in range(len(ubl)): d.append(np.append(ubl[i],0))
        d = np.array(d)
        m1 = -a.dot(la.pinv(np.transpose(a).dot(a))).dot(np.transpose(a))
        m2 = d.dot(la.pinv(np.transpose(a).dot(a))).dot(np.transpose(a))
        degenM = np.append(m1,m2,axis=0)
        #####################################################################################
        #A: A matrix for logcal amplitude
        A = np.zeros([len(crosspair),nAntenna+len(ubl)])
        for i in range(len(crosspair)):
            A[i][crosspair[i][0]] = 1
            A[i][crosspair[i][1]] = 1
            A[i][nAntenna+bltoubl[i]] = 1
        A = sps.csr_matrix(A)
        #################################################################################
        #B: B matrix for logcal phase
        B = np.zeros([len(crosspair),nAntenna+len(ubl)])
        for i in range(len(crosspair)):
            B[i][crosspair[i][0]] = reverse[i]*1 # XXX wha?
            B[i][crosspair[i][1]] = reverse[i]*-1 # XXX wha?
            B[i][nAntenna+bltoubl[i]] = 1
        B = sps.csr_matrix(B)
        ############################################################################
        #create info dictionary
        i = RedundantInfo(threshold=self.threshold)
        i.nAntenna = nAntenna
        i.nUBL = nUBL
        i.nBaseline = nBaseline
        i.subsetant = subsetant
        i.antloc = antloc
        i.subsetbl = subsetbl
        i.ubl = ubl
        i.bltoubl = bltoubl
        i.reversed = reverse
        i.reversedauto = reversedauto
        i.autoindex = autoindex
        i.crossindex = crossindex
        #i.ncross = ncross
        i.bl2d = bl2d
        i.ublcount = ublcount
        i.ublindex = ublindex
        i.bl1dmatrix = bl1dmatrix
        i.degenM = degenM
        i.A = A
        i.B = B
        i.totalVisibilityId = self.totalVisibilityId
        i.update()
        return i
    
