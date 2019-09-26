# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 15:34:34 2016

@author: Lina492375qW1188
"""

ATOM_NUMBER = 60150

def ReadAngle(filename):
    f = open('.\%s' %filename, 'r')
    adata = []
    for line in f:
        if 'ITEM: ATOMS id c_3 c_4 ' in line: 
            for line in f:
                adata.append(tuple([float(cell) for cell in line.split()]))
    f.close()
    adata = sorted(adata, key=lambda adata: adata[0]) # sort=排序 data
    return adata

def ReadUnrolled(filename):
    f = open('.\%s' %filename, 'r')
    unrolled=[]
    for line in f:
        if 'Atoms' in line:
            for i, line in enumerate(f): #enumerate = 照順序列出來
                unrolled.append(tuple([float(cell) for cell in line.split()]))
                if i > ATOM_NUMBER:
                    break
    f.close()

    del unrolled[0]
    return unrolled



def flatten(unrolled, data_S1, data_S, xmin, xmax):
    # unrolled: row1[3] = x coordinate
    #           row1[4] = y coordinate
    # data_S1:  row2[1] = eangle/atom at step 1
    #           row2[2] = ebond/atom at step 1
    # data_S:   row3[1] = eangle/atom at step S
    # data_S:   row3[2] = ebond/atom at step S
    separation = ([tuple([row1[3], row1[4], row2[1], row2[2], row3[1], row3[2]])
                     for (row1, row2, row3) in zip(unrolled, data_S1, data_S) 
                     if row1[3] > xmin and row1[3] < xmax])
    
    x = [row[0] for row in separation]
    y = [row[1] for row in separation]
    angle = [row[4]-row[2] for row in separation]
    bond = [row[5]-row[3] for row in separation]

    flatten_coord = list(zip(x, y, angle, bond))
    return flatten_coord

def OnlyCreases(angle_emap, bond_emap, threshold, edge, edge_threshold):
    
    L_a = len(angle_emap)
    
    for index, row_a, row_b in zip(range(L_a), angle_emap, bond_emap):
        
        if row_a[2] < threshold:
            angle_emap[index] = (row_a[0], row_a[1], 0.0)
            bond_emap[index] = (row_b[0], row_b[1], 0.0)
            
        elif abs(row_a[0]) > edge and row_a[2] < edge_threshold:
                angle_emap[index] = (row_a[0], row_a[1], 0.0)
                bond_emap[index] = (row_b[0], row_b[1], 0.0)
                
    return angle_emap, bond_emap

def RemoveRegion(angle_emap, bond_emap, xmin, xmax, ymin, ymax):
    
    L_a = len(angle_emap)
    
    for index, row_a, row_b in zip(range(L_a), angle_emap, bond_emap):
        
        if (row_a[0] < xmax and row_a[0] > xmin and 
            row_a[1] < ymax and row_a[1] > ymin):
            angle_emap[index] = (row_a[0], row_a[1], 0.0)
            bond_emap[index] = (row_b[0], row_b[1], 0.0)
            
    return angle_emap, bond_emap

# Data Input
output_S1 = ReadAngle('pe0.L360_K20000_A75_T0_V0.1_mkT_cy_Tn0_S1_a_r0(0)') #step0
output_S = ReadAngle('pe0.L190_K20000_A75_T1_V0.1_mkT_cy_Tn_S76670_a_r0(8)Es') #analysis
unrolled = ReadUnrolled('unrolled60150')

# Set parameters
Lside = 75 # left side
Rside = -75 # right side
darkness = 0.2
brightness = 0.5
magnification = 0.4
Nx = 100
Ny = 150

# Seperate angle and bond energy map.
x, y, angle, bond = zip(*flatten(unrolled, output_S1, output_S, Rside, Lside))
angle_tot = list(zip(x, y, angle))
bond_tot = list(zip(x, y, bond))

import PlotCollection
intensity = PlotCollection.Intensity2D(darkness, brightness, magnification)

# Plot total angle and bond energy map.
#intensity.ScipyBinned(angle_tot, Nx, Ny) 
#intensity.ScipyBinned(bond_tot , Nx, Ny) 

#eangle_tot = sum([row[2] for row in flat_angle_tot]) # calculate total angle energy.
#ebond_tot = sum([row[2] for row in flat_bond_tot])   # calculate total bond energy.



# Creases Angle and Bond Energy:
# OnlyCreases would cut those useless part first.
# useless part: including those ones with too high energy or where you don't want.
# Remember that you need this command to create angle_crease and bond_crease in
# order to use all the following functions.

angle_crease, bond_crease = OnlyCreases(angle_tot, bond_tot, 0, 100, 0)

#intensity.ScipyBinned(angle_crease, Nx, Ny) 
#intensity.ScipyBinned(bond_crease, Nx, Ny)

#eangle_crease = sum([row[2] for row in angle_crease])
#ebond_crease = sum([row[2] for row in bond_crease])

# return the data of Energy map
#bs = intensity.EnergyMap(angle_crease, Nx, Ny)


###############################################################################
# separate x, y, angle energy and bond energy from angle_crease and bond_crease.
# 
x, y, a_crease = zip(*angle_crease)
x, y, b_crease = zip(*bond_crease)
ab = [a+b for a, b in zip(a_crease, b_crease)]
ab_crease = list(zip(x, y, ab))
bs = intensity.EnergyMap(ab_crease, Nx, Ny)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

bz = bs.statistic

norm = (cm.colors.Normalize(vmin=np.nanmin(bz)*darkness,
                            vmax=np.nanmax(bz)*brightness))

fig = plt.figure()

(plt.imshow(bz, extent=(np.amin(bs.y_edge), np.amax(bs.y_edge),
                                np.amin(bs.x_edge), np.amax(bs.x_edge)),
                                interpolation='sinc', origin='lower',
                                aspect='auto', cmap=cm.afmhot, norm=norm))
#, origin='lower'
d = (np.nanmax(bz)-np.nanmin(bz))/4
(plt.colorbar(ticks=[np.nanmin(bz), np.nanmin(bz)+d, np.nanmin(bz)+2*d, 
                     np.nanmin(bz)+3*d, np.nanmax(bz)]))

plt.xlabel('y')
plt.ylabel('x')

# FIND_CLOSEST_POINT
def Find(x, arr):
    """Find the point in arr which is closest to value."""
    ref = np.array([x]*len(arr))
    d = [(index, abs(diff)) for index, diff in enumerate(ref-arr)]    
    # sort array d and make the smallest difference be the first one. And then
    # return the index of it.
    i_smallest = sorted(d, key=lambda b:b[1])[0][0]
    return i_smallest

# DEFINE_POLYNOMIAL
def PolyFunc(press_arr):
    """press_arr should contain at least 2 elements. This function will return
       the polynomial function that connect all points in press_arr."""
    polyfunc_list = []
    # Reset polyfunc_list everytime for each call, which is used to avoid 
    # ambiguity and make each call of PolyFunc return list of polynomial 
    # functions that corresponds to input array press_arr.
    
    for i in range(0, len(press_arr)-1):
        x1, y1 = press_arr[i][2], press_arr[i][3]
        x2, y2 = press_arr[i+1][2], press_arr[i+1][3]
               
        if y1 != y2:
            m = (x2-x1)/(y2-y1) # slop
            b = x1 - m*y1       # intercept
            polyfunc_list.append(np.poly1d([m, b]))
        else:
            b = y1
            polyfunc_list.append(np.poly1d([b]))
        
    return polyfunc_list

def pos_on_pline(ends, polyfunc):
    #ends = [(x_i, y_i, xdata_i, ydata_i), 
    #        (x_{i+1}, y_{i+1}, xdata_{i+1}, ydata_{i+1})]
# Pixel:
    #x_begin, y_begin = ends[0][0], ends[0][1]
    #x_end, y_end = ends[1][0], ends[1][1]
    
    xdata_begin, ydata_begin = ends[0][2], ends[0][3]
    xdata_end, ydata_end = ends[1][2], ends[1][3]

# how many value in bs.x_edge between two point.
    n = abs(Find(ydata_end, bs.x_edge)-Find(ydata_begin, bs.x_edge))
# how many pixel between two points (b.t.p.).
    #n = int(abs(y_end - y_begin)) 
    dx = (xdata_end - xdata_begin)/n # the data's real difference b.t.p.
    dy = (ydata_end - ydata_begin)/n
    
    pos_list = [] # position list stores all the positions of point on line.
    for i in range(n):
        x_pos = polyfunc(ydata_begin + i*dy)  # xpos is belong to y_edge
        y_pos = ydata_begin + i*dy            # ypos is belong to x_edge
        #z_pos = bz[Find(y_pos, bs.x_edge)][Find(x_pos, bs.y_edge)]
        
        x_real = bs.x_edge[Find(y_pos, bs.x_edge)]
        y_real = bs.y_edge[Find(x_pos, bs.y_edge)]
        # x_real is in x_edge and y_real is in y_edge.
        pos_list.append((x_real, y_real))
        
    return pos_list


### GAUSSIAN_SUM
from scipy.optimize import curve_fit

def Gaussian(x, a, x0, w):
    e0 = -(x-x0)**2/(2*w**2)
    return a*np.exp(e0)

def Area(xc, yc):
    """xc is in x_edge and yc is in y_edge."""
    t = 0.5 # threshold of yc.

    # Y-fit:
    Y = bz[Find(xc, bs.x_edge)]
    dY = np.gradient(np.log10(Y))     # 1st derivative
    ddY = np.gradient(dY)             # 2nd derivative
    k = abs(ddY[Find(yc, bs.y_edge)]) # curvature at point nearest to yc
    w = np.sqrt(1/(2*k))              # the estimate of width
    
    lower = [0, yc-t, 0]
    higher = [np.inf, yc+t, w]
    
    height = bz[Find(xc, bs.x_edge)]
    
    ydiff = bs.y_edge[1]-bs.y_edge[0]
    horiz = np.delete(bs.y_edge, -1) + ydiff
    
    para = curve_fit(Gaussian, horiz, height, bounds=(lower, higher))
    a, x0, w = para[0]

    # X-fit:
    X = bz[:, Find(yc, bs.y_edge)]
    dX = np.gradient(np.log10(X))
    ddX = np.gradient(dX)
    kX = abs(ddX[Find(xc, bs.x_edge)])
    Xw = np.sqrt(1/(2*kX))
    
    Xlower = [a-t, xc-t, 0]
    Xhigher = [a+t, xc+t, Xw]
    
    Xheight = bz[:, Find(yc, bs.y_edge)]
    
    xdiff = bs.x_edge[1]-bs.x_edge[0]
    Xhoriz = np.delete(bs.x_edge, -1) + xdiff
    
    para = curve_fit(Gaussian, Xhoriz, Xheight, bounds=(Xlower, Xhigher))
    a, Xx0, Xw = para[0]

    #return np.sqrt(2*np.pi)*a*w
    return np.pi*Xw*w*a

### INTERACTIVE_PRESS
L = []
area = []
press_arr = []

def on_press(event, polyfunc = []):
    
    press = event.x, event.y, event.xdata, event.ydata
    press_arr.append(press)
    
    try:
        polyfunc = PolyFunc(press_arr)
        ends = [press_arr[-2], press_arr[-1]]
        
# Calculate length:
        Li = (np.sqrt((press_arr[-1][2]-press_arr[-2][2])**2+
                     (press_arr[-1][3]-press_arr[-2][3])**2))
        L.append(Li)
        
# Extract and draw points:
        pos_list = pos_on_pline(ends, polyfunc[-1])
        pos_x, pos_y = zip(*pos_list)
        plt.plot(pos_y, pos_x, marker = 'p', markersize = 6, color = 'b')
        plt.xlim((np.amin(bs.y_edge), np.amax(bs.y_edge)))
        plt.ylim((np.amin(bs.x_edge), np.amax(bs.x_edge)))
        fig.canvas.draw()

# Calculate area:        
        area.append(sum([Area(xc, yc) for xc, yc in pos_list]))

    except:
        print("Only 1 point in press_arr!")

def on_key(event): # 空白鍵會執行兩次... Don't ask me why!!!

    print("=================================================")
    print("The length of all line segment:", sum(L))
    print("The sum of all line segment:", sum(area))
    print("=================================================")

    f = open('.\data.txt','a')
    f.write(str(sum(L))+' '+str(sum(area))+'\n')
    f.close
    
    del L[:]         # Reset L
    del area[:]      # Reset area
    del press_arr[:] # Reset press_arr


cidpress = fig.canvas.mpl_connect('button_press_event', on_press)
cidkey = fig.canvas.mpl_connect('key_press_event', on_key)

plt.show()




############################################################################### 
#CHECK_w

option = 'a'

if option == 'b':
    import numpy as np
    import matplotlib.pyplot as plt
    
    xc = 53.73 # the position in the axial direction of paper roll.
    yc = 5.0
    
# Use curvature to find the estimate of w:
    Y = bz[Find(xc, bs.x_edge)]
    dY = np.gradient(np.log10(Y))
    ddY = np.gradient(dY)
    k = abs(ddY[Find(yc, bs.y_edge)])
    w = np.sqrt(1/(2*k))
    
    lower = [0, yc-0.5, 0]
    higher = [np.inf, yc+0.5, w+0.5]
    
    ydiff = bs.y_edge[1]-bs.y_edge[0]
    X = bs.y_edge[0:300] + ydiff
    
    p0 = [1, yc, 1]
    para = curve_fit(Gaussian, X, Y, p0=p0, bounds=(lower, higher))
    a, x0, w = para[0]
    Y_fit = np.array([Gaussian(Xi, a, x0, w) for Xi in X])
    
    plt.plot(X, Y)
    plt.plot(X, Y_fit)
    plt.show()
