# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 13:07:15 2016

@author: 510a
"""

import math
import numpy as np

N=2
K=1
n=[0,1,2]

def fun(list_n):
    coord=[]
    for i in list_n:
        coord.append(tuple([math.cos(math.pi*i/N), math.sin(math.pi*i/N)]))
        
    return coord

coord = fun(n)
coord = np.array(coord)

v1 = coord[0]-coord[1]
v2 = coord[2]-coord[1]

nv1 = np.linalg.norm(v1)
nv2 = np.linalg.norm(v2)

cos_angle = np.dot(v1,v2)/(nv1*nv2)

eangle=K*(1-cos_angle)*N
