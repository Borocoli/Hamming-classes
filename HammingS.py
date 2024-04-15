import numpy as np

class HammingS():
    def __init__(self, n, k):
        self.n = n
        self.k = k
        self.m = n-k
        l = []
        im = []
        q = []
        for i in range(1, 2**self.m):
            v = to_bin(i, self.m)
            l.append(v) 
            if np.sum(v) == 1:
                im = [v] + im
            else:
                q.append(v)
        self.H = np.concatenate(l).reshape(2**self.m - 1, self.m).T
        I = np.concatenate(im).reshape(self.m, self.m).T
        Q = np.concatenate(q).reshape(2**self.m - 1 - self.m, self.m).T
        self.H = np.concatenate((I, Q), axis=1)
        self.G = np.concatenate((Q.T, np.identity(self.k, dtype=int)), axis=1)

    def encode(self,el, canonic=True):
        if type(el) == str or type(el) == list:
            el = [int(c) for c in el]
        i = np.asarray(el, dtype=int)
        if len(i) < self.k:
            i = np.concatenate((i, np.zeros(self.k - len(i), dtype=int)), dtype=int)
        if canonic:
            return i@self.G%2
        else:
            vi = i@self.G%2
            v = [None]*(2**self.m-1)
            for a in range(0, self.m):
                v[2**a -1] = vi[self.m -1 - a]
            k = self.m
            for a in range(len(v)):
                if v[a] == None:
                    v[a] = vi[k]
                    k+= 1
            return v

    def decode(self, el, canonic=True):
        if type(el) == str or type(el) == list:
            el = [int(c) for c in el]
        v = np.asarray(el, dtype=int)
        if len(v) < self.n:
            v = np.concatenate((v, np.zeros(self.n - len(v), dtype=int)), dtype=int)
        if not canonic:
            c = []
            j = []
            for i in range(len(v)):
                 if np.sum(to_bin(i+1)) == 1:
                     c = [v[i]] + c
                 else:
                     j.append(v[i])
            v = c + j
            v = np.array(v, dtype=int)
        z = self.H@v.T%2
        if np.sum(z):
            Ht = self.H.T
            for i in range(len(Ht)):
                if np.all(z == Ht[i]):
                    v[i] = (v[i] + 1)%2
        rez = np.asarray(v[self.m:])
        return rez%2
                                 

''
