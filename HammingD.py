import numpy as np
from sympy import *

class HammingD():
    def __init__(self, n, k):
        self.n = n
        self.k = k
        self.m = n-k-1
        l = []
        im = []
        q = []
        for i in range(1, 2**self.m):
            v = to_bin(i, self.m)
            l.append(v) 
        self.H = np.concatenate(l).reshape(2**self.m - 1, self.m).T
        self.H = np.concatenate((self.H, np.zeros((self.m, 1), dtype=int)), axis=1)
        self.H = np.concatenate((self.H, np.ones((1, 2**self.m), dtype=int)))
        self.m += 1
        v = []
        c = []
        j = []
        for i in range(self.n):
            a = ''
            if np.sum(to_bin(i+1)) == 1:
                a = Symbol(f'c{1 + int(np.log2(i+1))}')
                c.append(a)
            else:
                a = Symbol(f'i{i - int(np.log2(i+1))}')
                j.append(a)
            v.append(a)
        self.j = j
        V = Matrix(v)
        H = Matrix(self.H)
        jj = [1, 0, 1, 1]
        A = H*V
        T = (A).subs([(j[i], jj[i]) for i in range(self.k)])
        z = [] 
        for i in A:
            i = i.subs(z)
            unknown = set(c)&i.free_symbols
            T = solve(i, [e for e in unknown])
            V = V.subs([(list(unknown)[k], T[k]) for k in range(len(T))])
            A = A.subs([(list(unknown)[k], T[k]) for k in range(len(T))])
            z +=[(list(unknown)[k], T[k]) for k in range(len(T))] 
            tmp = [list(unknown)[k] for k in range(len(T))]
            for k in tmp:
                c.remove(k)
        self.V = V


    def encode(self,el, canonic=True):
        if type(el) == str or type(el) == list:
            el = [int(c) for c in el]
        i = np.asarray(el, dtype=int)
        if len(i) < self.k:
            i = np.concatenate((i, np.zeros(self.k - len(i), dtype=int)), dtype=int)
        T = [(self.j[k], i[k]) for k in range(self.k)]
        v = self.V.subs(T)%2
        v = np.asarray(v.T, dtype=int).reshape(len(v))
        return v

    def decode(self, el, canonic=True):
        if type(el) == str or type(el) == list:
            el = [int(c) for c in el]
        v = np.asarray(el, dtype=int)
        if len(v) < self.n:
            v = np.concatenate((v, np.zeros(self.n - len(v), dtype=int)), dtype=int)

        v = Matrix(v)
        z = self.H*v%2
        z = np.asarray(z).reshape(len(z))
        if np.any(z):
            Ht = self.H.T
            found = False
            for i in range(len(Ht)):
                if np.all(z == Ht[i]):
                    v[i] = (v[i] + 1)%2
                    found = True
                    break
            if not found:
                return -1
        v = self.V - v
        sol = []
        for line in v:
            unknown = set(self.j)&line.free_symbols
            if len(unknown) != 1:
                continue
            sol.append(solve(line, list(unknown))[0]%2)
        sol = np.asarray(sol, dtype=int)
        return sol


