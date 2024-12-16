#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import timeit

my_setup = '''
import numpy as np
from generalizedplanck.components import (GeneralizedPlanck, Reflectance, 
                                          IdealSqrtAbsorption, UrbachTail)
from generalizedplanck.utils import get_nk
from generalizedplanck.data import gaas_cl

s = gaas_cl()

n1 = get_nk(n=1)
n2 = get_nk(n=3)

r = Reflectance(theta=0, n1=n1, n2=n2, pol=None) 
ideal_abs = IdealSqrtAbsorption(Eg=1.42, E0=1.6, a0=14800)
tail = UrbachTail(g=0.015)

genp_anal = GeneralizedPlanck(Eg=1.42, g=0.01, p=0.8, T=300, d=150, 
                              Efv=0, Efc=0.1,
                              reflectance=r,
                              ideal_abs_coeff=ideal_abs,
                              tail=tail,
                              analytical=True)

genp_brute = GeneralizedPlanck(Eg=1.42, g=0.01, p=0.8, T=300, d=150, 
                              Efv=0, Efc=0.1,
                              reflectance=r,
                              ideal_abs_coeff=ideal_abs,
                              tail=tail,
                              analytical=False)

en = np.linspace(1.3, 1.6, 1000)
'''

brutal_conv = '''genp_brute.abs_coeff_tail(en)'''
analytical_conv = '''genp_anal.abs_coeff_tail(en)'''

analytical_time = timeit.timeit(setup=my_setup, stmt=analytical_conv, 
                                number=10000)/10000
brutal_time = timeit.timeit(setup=my_setup, stmt=brutal_conv, number=10)/10

print(f'Analytical time : {analytical_time} s')
print(f'Brute convolution time : {brutal_time} s')