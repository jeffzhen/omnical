import numpy as np
import numpy.linalg as la
import scipy.sparse as sps
from info import RedundantInfo

class ArrayInfo:
    '''Store information about an antenna array needed for computing redundancy and indexing matrices.'''
    def __init__(self, nTotalAnt, badAntenna=[], badUBLpair=[]):
        self.nTotalAnt = nTotalAnt
        self.nTotalBaselineAuto = (nTotalAnt + 1) * nTotalAnt / 2
        self.nTotalBaselineCross = (nTotalAnt - 1) * nTotalAnt / 2
        self.antennaLocation = np.zeros((nTotalAnt, 3))
        side = int(nTotalAnt**.5)
        for a in range(nTotalAnt): self.antennaLocation[a] = np.array([a/side, a%side, 0])
        self.antennaLocationTolerance = 1e-6
        self.badAntenna = badAntenna
        self.badUBLpair = badUBLpair
        #PAPER miriad convention by default
        self.totalVisibilityId = np.concatenate([[[i,j] for i in range(j+1)] for j in range(nTotalAnt)])
        self._gen_totalVisibilityId_dic()
        self.totalVisibilityUBL = None
    def _gen_totalVisibilityId_dic(self):
        self.totalVisibilityId_dic = {}
        for i, (a1,a2) in enumerate(self.totalVisibilityId): self.totalVisibilityId_dic[(a1,a2)] = i
    def get_baseline(self,bl):
        '''inverse function of totalVisibilityId, calculate the bl index from 
        the antenna pair. It allows flipping of a1 and a2, will return same result'''
        bl = tuple(bl)
        try: return self.totalVisibilityId_dic[bl]
        except(KeyError): pass
        try: return self.totalVisibilityId_dic[bl[::-1]]
        except(KeyError): return None
    def read_arrayinfo(self, arrayinfopath, verbose=False):
        '''array info is the minimum set of information to uniquely describe a 
        redundant array, and is needed to compute redundant info. It includes, 
        in each line, bad antenna indices, bad unique bl indices, tolerance 
        of error when checking redundancy, antenna locations, and visibility's 
        antenna pairing conventions. Unlike redundant info which is a self-contained 
        dictionary, items in array info each have their own fields in the instance.'''
        if verbose: print "Reading", arrayinfopath
        with open(arrayinfopath) as f: rawinfo = [[float(x) for x in line.split()] for line in f]
        self.badAntenna = np.array(rawinfo[0], dtype=np.int)
        if self.badAntenna[0] < 0: self.badAntenna = np.zeros(0) # XXX special significance for < 0?
        rawpair = np.array(rawinfo[1], dtype=np.int)
        if rawpair.shape[0] == 0 or rawpair.shape[0] % 2 != 0 or rawpair.min() < 0: # XXX shouldn't accept bad states
            self.badUBLpair = np.array([])
        else: self.badUBLpair = np.reshape(rawpair,(len(rawpair)/2,2))
        self.antennaLocationTolerance = rawinfo[2][0]
        for a in range(len(self.antennaLocation)):
            assert(len(rawinfo[a+3]) == 3)
            self.antennaLocation[a] = np.array(rawinfo[a+3])
        bl = 0
        vis_id = []
        max_bl_cnt = self.nTotalAnt * (self.nTotalAnt + 1) / 2
        maxline = len(rawinfo)
        while len(rawinfo[bl + 3 + len(self.antennaLocation)]) == 2: # XXX don't like while loop
            assert(bl < max_bl_cnt)
            vis_id.append(np.array(rawinfo[bl + 3 + len(self.antennaLocation)], dtype=np.int))
            bl += 1
            if bl + 3 + len(self.antennaLocation) >= maxline: break
        self.totalVisibilityId = np.array(vis_id, dtype=np.int)
        self._gen_totalVisibilityId_dic()
    def compute_UBL(self,tolerance = 0.1):
        '''XXX DOCSTRING'''
        if tolerance == 0:
            tolerance = np.min(np.linalg.norm(np.array(self.antennaLocation) - self.antennaLocation[0], axis=1)) / 1e6
        ubl = {}
        for bl, (a1, a2) in enumerate(self.totalVisibilityId):
            if a1 != a2 and a1 not in self.badAntenna and a2 not in self.badAntenna:
                loc_tuple = tuple(np.round((self.antennaLocation[a2] - self.antennaLocation[a1]) / float(tolerance)) * tolerance)
                neg_loc_tuple = tuple(np.round((self.antennaLocation[a1] - self.antennaLocation[a2]) / float(tolerance)) * tolerance)
                if loc_tuple in ubl: ubl[loc_tuple].add(bl + 1)
                elif neg_loc_tuple in ubl: ubl[neg_loc_tuple].add(- bl - 1)
                else:
                    if loc_tuple[0] >= 0: ubl[loc_tuple] = set([bl + 1])
                    else: ubl[neg_loc_tuple] = set([-bl - 1])
        #calculate actual average of the gridded bls vectors to get an accurate representation of the ubl vector
        ubl_vec = np.zeros((len(ubl), 3))
        self.totalVisibilityUBL = {}
        ublcount = np.zeros(len(ubl))
        for u, grid_ubl_vec in enumerate(ubl):
            for bl in ubl[grid_ubl_vec]:
                assert bl != 0
                a1, a2 = self.totalVisibilityId[abs(bl) - 1]
                if bl > 0: ubl_vec[u] = ubl_vec[u] + self.antennaLocation[a2] - self.antennaLocation[a1]
                else: ubl_vec[u] = ubl_vec[u] + self.antennaLocation[a1] - self.antennaLocation[a2]
                self.totalVisibilityUBL[(a1, a2)] = u
            ublcount[u] = len(ubl[grid_ubl_vec])
            ubl_vec[u] = ubl_vec[u] / ublcount[u]
        reorder = (ubl_vec[:,1]*1e9 + ubl_vec[:,0]).argsort()
        rereorder = reorder.argsort()
        for key in self.totalVisibilityUBL:
            self.totalVisibilityUBL[key] = rereorder[self.totalVisibilityUBL[key]]
        ubl_vec = ubl_vec[reorder]
        #now I need to deal with the fact that no matter how coarse my grid is, it's possible for a single group of ubl to fall into two adjacent grids. So I'm going to check if any of the final ubl vectors are seperated by less than tolerance. If so, merge them
        ublmap = {}
        for u1 in range(len(ubl_vec)):
            for u2 in range(u1):
                if la.norm(ubl_vec[u2] - ubl_vec[u1]) < tolerance or la.norm(ubl_vec[u2] + ubl_vec[u1]) < tolerance:
                    ublmap[u1] = u2
                    ubl_vec[u2] = (ubl_vec[u1] * ublcount[u1] + ubl_vec[u2] * ublcount[u2]) / (ublcount[u1] + ublcount[u2])
                    break
            ublmap[u1] = u1
        merged_ubl_vec = []
        for u in range(len(ubl_vec)):
            if ublmap[u] == u:
                merged_ubl_vec.append(ubl_vec[u])
                ublmap[u] = len(merged_ubl_vec) - 1
            else: ublmap[u] = ublmap[ublmap[u]]
        merged_ubl_vec = np.array(merged_ubl_vec)
        for key in self.totalVisibilityUBL:
            self.totalVisibilityUBL[key] = ublmap[self.totalVisibilityUBL[key]]
        return ubl_vec
    def compute_redundantinfo(self, arrayinfoPath=None, verbose=False, badAntenna=[], badUBLpair=[], antennaLocationTolerance=1e-6):
        '''Use provided antenna locations (in arrayinfoPath) to derive redundancy equations'''
        # XXX could these be set somewhere else so they aren't passed in?
        self.antennaLocationTolerance = tol = antennaLocationTolerance
        self.badAntenna += badAntenna
        self.badUBLpair += badUBLpair
        if arrayinfoPath is not None: self.read_arrayinfo(arrayinfoPath)
        info = RedundantInfo()
        # exclude bad antennas
        info['subsetant'] = subsetant = np.array([i for i in xrange(self.antennaLocation.shape[0]) 
                if i not in self.badAntenna], dtype=np.int32)
        info['nAntenna'] = nAntenna = len(subsetant) # XXX maybe have C api automatically infer this
        info['antloc'] = antloc = np.array([self.antennaLocation[i] for i in subsetant], dtype=np.float32)
        ublall = self.compute_UBL(tol)
        #delete the bad ubl's
        badUBL = {}
        def dis(a1,a2): return np.linalg.norm(a1-a2)
        for a1,a2 in self.badUBLpair:
            bl = self.antennaLocation[a1] - self.antennaLocation[a2]
            for i,ubl in enumerate(ublall):
                if dis(bl,ubl) < tol or dis(bl,-ubl) < tol: badUBL[i] = None
        ubl2goodubl = {}
        def f(i,u):
            ubl2goodubl[i] = len(ubl2goodubl)
            return u
        info['ubl'] = ubl = np.array([f(i,u) for i,u in enumerate(ublall) if not badUBL.has_key(i)], dtype=np.float32)
        for k in badUBL: ubl2goodubl[k] = -1
        #info['nUBL'] = nUBL = ubl.shape[0] # XXX maybe have C api automatically infer this
        nUBL = ubl.shape[0] # XXX maybe have C api automatically infer this
        badubl = [ublall[i] for i in badUBL]
        #find nBaseline (include auto bls) and subsetbl
        #bl2d:  from 1d bl index to a pair of antenna numbers
        bl2d = [] # XXX cleaner way to do this?
        for i,ai in enumerate(antloc):
            for j,aj in enumerate(antloc[:i+1]):
                blij = ai - aj
                flag = False
                for bl in badubl:
                    if dis(blij,bl) < tol or dis(blij,-bl) < tol:
                        flag = True
                        break
                if not flag: bl2d.append((i,j))
        # exclude pairs that are not in totalVisibilityId
        tmp = []
        for p in bl2d:
            bl = (subsetant[p[0]],subsetant[p[1]])
            if self.totalVisibilityId_dic.has_key(bl): tmp.append(p)
            elif self.totalVisibilityId_dic.has_key(bl[::-1]): tmp.append(p[::-1])
        #info['bl2d'] = bl2d = np.array(tmp, dtype=np.int32)
        bl2d = np.array(tmp, dtype=np.int32)
        crossindex = np.array([i for i,p in enumerate(bl2d) if p[0] != p[1]], dtype=np.int32)
        nBaseline = len(bl2d)
        bl2d = bl2d[crossindex] # make bl2d only hold crosscorrelations
        info['nBaseline'] = len(bl2d) # XXX maybe have C api infer this
        # from a pair of good antenna index to bl index
        info['subsetbl'] = np.array([self.get_baseline([subsetant[bl[0]],subsetant[bl[1]]]) 
                for bl in bl2d], dtype=np.int32)
        #bltoubl: cross bl number to ubl index
        def findublindex(p1,p2):
            a1,a2 = subsetant[p1],subsetant[p2]
            if (a1,a2) in self.totalVisibilityUBL: return ubl2goodubl[self.totalVisibilityUBL[(a1,a2)]]
        info['bltoubl'] = bltoubl = np.array([findublindex(*p) for p in bl2d if p[0] != p[1]], dtype=np.int32)
        #reversed:   cross only bl if reversed -1, otherwise 1
        crosspair = [p for p in bl2d if p[0] != p[1]]
        reverse = []
        for k,cpk in enumerate(crosspair):
            bl = antloc[cpk[0]] - antloc[cpk[1]]
            if dis(bl,ubl[bltoubl[k]]) < tol: reverse.append(-1)
            elif dis(bl,-ubl[bltoubl[k]]) < tol: reverse.append(1)
            else : raise ValueError('bltoubl[%d] points to wrong ubl index' % (k))
        reverse = np.array(reverse, dtype=np.int32)
        #info['reversed'] = np.ones_like(reverse)
        info._reversed = reverse # XXX store this to remember what we did
        bl2d0 = np.where(reverse == 1, bl2d[:,0], bl2d[:,1])
        bl2d1 = np.where(reverse == 1, bl2d[:,1], bl2d[:,0])
        bl2d[:,0],bl2d[:,1] = bl2d0,bl2d1
        crosspair = [p for p in bl2d if p[0] != p[1]] # recompute crosspair for reversed indices
        info.bl2d = bl2d
        #reversedauto = np.arange(info['nBaseline'], dtype=np.int32)
        #for i in autoindex: reversedauto[i] = 1
        #ublcount:  for each ubl, the number of good cross bls corresponding to it
        cnt = {}
        for bl in bltoubl: cnt[bl] = cnt.get(bl,0) + 1
        info['ublcount'] = np.array([cnt[i] for i in range(nUBL)], dtype=np.int32)
        #ublindex:  //for each ubl, the vector<int> contains (ant1, ant2, crossbl)
        cnt = {}
        for i,(a1,a2) in enumerate(crosspair): cnt[bltoubl[i]] = cnt.get(bltoubl[i],[]) + [[a1,a2,i]]
        #info['ublindex'] = ublindex = np.concatenate([np.array(cnt[i],dtype=np.int32) for i in range(nUBL)])
        ublindex = np.concatenate([np.array(cnt[i],dtype=np.int32) for i in range(nUBL)])
        newind = np.arange(nBaseline)[crossindex] = np.arange(crossindex.size, dtype=np.int32)
        ublindex[2] = newind[ublindex[2]]
        info.ublindex = ublindex
        #bl1dmatrix: a symmetric matrix where col/row numbers index ants and entries are bl index (no auto corr)
        # XXX don't like 2**31-1.  whence this number?
        bl1dmatrix = (2**31-1) * np.ones((nAntenna,nAntenna),dtype=np.int32)
        for i,cp in enumerate(crosspair): bl1dmatrix[cp[1],cp[0]], bl1dmatrix[cp[0],cp[1]] = i,i
        info['bl1dmatrix'] = bl1dmatrix
        #degenM:
        a = np.array([np.append(ai,1) for ai in antloc], dtype=np.float32)
        d = np.array([np.append(ubli,0) for ubli in ubl], dtype=np.float32)
        m1 = -a.dot(la.pinv(a.T.dot(a))).dot(a.T)
        m2 = d.dot(la.pinv(a.T.dot(a))).dot(a.T)
        info['degenM'] = np.append(m1,m2,axis=0)
        #A: A matrix for logcal amplitude
        A = np.zeros((len(crosspair),nAntenna+nUBL))
        for i,cp in enumerate(crosspair): A[i,cp[0]], A[i,cp[1]], A[i,nAntenna+bltoubl[i]] = 1,1,1
        info['At'] = sps.csr_matrix(A).T
        #B: B matrix for logcal phase
        B = np.zeros((len(crosspair),nAntenna+nUBL))
        #for i,cp in enumerate(crosspair): B[i,cp[0]], B[i,cp[1]], B[i,nAntenna+bltoubl[i]] = -reverse[i],reverse[i],1
        for i,cp in enumerate(crosspair): B[i,cp[0]], B[i,cp[1]], B[i,nAntenna+bltoubl[i]] = -1,1,1
        info['Bt'] = sps.csr_matrix(B).T
        info.update()
        return info
