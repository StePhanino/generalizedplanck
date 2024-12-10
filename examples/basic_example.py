#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 15:41:31 2024

@author: stefano
"""

import hyperspy.api as hs

from generalizedplanck.components import (GeneralizedPlanck, Reflectance, 
                                          IdealSqrtAbsorption, UrbachTail)
from generalizedplanck.utils import get_nk
from generalizedplanck.data import gaas_cl

s = gaas_cl()

n1 = get_nk(n=1)
n2 = get_nk(shelf='main', book='GaAs', page='Papatryfonos')

gaas_r = Reflectance(theta=0, n1=n1, n2=n2, pol=None) 
gaas_ideal_abs = IdealSqrtAbsorption(Eg=1.42, E0=1.6, a0=14800)
gaas_tail = UrbachTail(g=0.015)

gaas_genp = GeneralizedPlanck(Eg=1.42, g=0.01, p=0.8, T=300, d=150, 
                              Efv=0, Efc=0.1,
                              reflectance=gaas_r,
                              ideal_abs_coeff=gaas_ideal_abs,
                              tail=gaas_tail)


