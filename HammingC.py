import numpy as np
from sympy import *

class HammingC():
    def __init__(self, n, i,  poly):
        self.g = parse_expr(poly, evaluate=True).as_poly(domain=GF(2))
        self.i = i
        self.c = n-i
        self.n = n
        
        x = symbols('x')
        self.h = parse_expr(f'x**{n} + 1', evaluate=True).as_poly()
        self.h, r = div(self.h,self.g)
        rez = self.h.all_coeffs()
        rez = np.array([intb(i) for i in rez], dtype=object) + 0
        rez = np.polyadd(np.zeros(n, dtype=int), np.flip(rez))
        self.H = np.zeros((self.c, n), dtype=int)
        for i in range(self.c):
            self.H[i] = rez
            rez = np.roll(rez, -1)
    
    def encode(self, cod):
        x = symbols('x')
        j = Poly(cod, x)
        mare = parse_expr(f'x**{self.c}', evaluate=True).as_poly()
        j = j*mare
        c = j%self.g
        rez = j + c
        rez = rez.all_coeffs()
        rez = np.array([intb(i) for i in rez], dtype=object) + 0
        return rez

    def decode(self, cod):
        x = symbols('x')
        p = Poly(cod, x)
        
        z = p%self.g
        zz = z.all_coeffs()
        #print(zz)
        #print(z - Poly([0], x))
        if not (z - Poly([0], x)) == Poly([0], x):
            a = parse_expr(f'x', evaluate=True).as_poly(domain=GF(2))
            print(a)
            for i in range(self.n):
                k = parse_expr(f'x**{i+1}', evaluate=True).as_poly(domain=GF(2))
                k = (k/a)%self.g
                print(k)

                if (z - p) == Poly([0], x):
                    print(i)
                    cod[i] = (cod[i]+1)%2
        return cod[:self.i]

