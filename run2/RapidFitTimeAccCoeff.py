import sys
import numpy as np 
import pylab as pl
import numpy.random as rand  # Random number generation module

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json

# Make it work for Python 2+3 and with Unicode
import io
try:
    to_unicode = unicode
except NameError:
    to_unicode = str

# Read JSON file
with open('/home/s1668588/Bs2JpsiPhi-Run2/ANALYSIS/analysis/fitinputs/fitinputs_st015_trigCat.json') as data_file:
    data_loaded = json.load(data_file)
    data_file.close() # to enable writing to later on

# These are the P,Q,R,S terms
def P(i):
    d1 = knots[i+1]-knots[i-2]
    d2 = knots[i+1]-knots[i-1]
    d3 = knots[i+1]-knots[i]
    #return 1.
    return d1*d2*d3

def Q(i):
    d1 = knots[i+2]-knots[i-1]
    d2 = knots[i+1]-knots[i-1]
    d3 = knots[i+1]-knots[i]
    #return 1.
    return d1*d2*d3

def R(i):
    d1 = knots[i+2]-knots[i]
    d2 = knots[i+2]-knots[i-1]
    d3 = knots[i+1]-knots[i]
    #return 1.
    return d1*d2*d3

def S(i):
    d1 = knots[i+3]-knots[i]
    d2 = knots[i+2]-knots[i]
    d3 = knots[i+1]-knots[i]
    #return 1.
    return d1*d2*d3

# This calculates the u^n coefficients pn for (u-u1)(u-u2)(u-u3)
def powerCoeffs( i1, i2, i3):
    p0=-(knots[i1]*knots[i2]*knots[i3])
    p1=(knots[i1]*knots[i2]+knots[i2]*knots[i3]+knots[i3]*knots[i1])
    p2=-(knots[i1]+knots[i2]+knots[i3])
    p3=1.
    return np.array([p0,p1,p2,p3])

# This calculates the value of (u-u1)(u-u2)(u-u3)
def curve( u, i1, i2, i3):
    return (u-knots[i1])*(u-knots[i2])*(u-knots[i3])

# Power coeffs and value of the "A" - "D" terms
def ApowerCoeffs(i):
    t1 = - powerCoeffs( i+1, i+1, i+1 ) / P(i)
    return t1
def Acurve(u, i):
    t1 = - curve(u, i+1, i+1, i+1 ) / P(i)
    return t1

def BpowerCoeffs(i):
    t1 = powerCoeffs( i-2, i+1, i+1 ) / P(i)
    t2 = powerCoeffs( i-1, i+1, i+2 ) / Q(i)
    t3 = powerCoeffs( i, i+2, i+2 ) / R(i)
    return t1+t2+t3
def Bcurve(u, i):
    t1 = curve( u, i-2, i+1, i+1 ) / P(i)
    t2 = curve( u, i-1, i+1, i+2 ) / Q(i)
    t3 = curve( u, i, i+2, i+2 ) / R(i)
    return t1+t2+t3

def CpowerCoeffs(i):
    t1 = - powerCoeffs( i-1, i-1, i+1 ) / Q(i)
    t2 = - powerCoeffs( i-1, i, i+2 ) / R(i)
    t3 = - powerCoeffs( i, i, i+3 ) / S(i)
    return t1+t2+t3
def Ccurve(u, i):
    t1 = - curve( u, i-1, i-1, i+1 ) / Q(i)
    t2 = - curve( u, i-1, i, i+2 ) / R(i)
    t3 = - curve( u, i, i, i+3 ) / S(i)
    return t1+t2+t3

def DpowerCoeffs(i):
    t1 = powerCoeffs( i, i, i ) / S(i)
    return t1
def Dcurve(u, i):
    t1 = curve( u, i, i, i ) / S(i)
    return t1

# This creates the sum of the power coeffs for knot interval j
def SumpowerCoeffs(bcoeffs,j):
    Alist = ApowerCoeffs( j+istart )
    Blist = BpowerCoeffs( j+istart )
    Clist = CpowerCoeffs( j+istart )
    Dlist = DpowerCoeffs( j+istart )
    Sumlist = bcoeffs[j]*Alist + bcoeffs[j+1]*Blist + bcoeffs[j+2]*Clist + bcoeffs[j+3]*Dlist
    return Sumlist

# This returns the value od S(u) for a given u in knot interval j
def SumCurve(bcoeffs, u, j):
    Ac = Acurve( u, j+istart )
    Bc = Bcurve( u, j+istart )
    Cc = Ccurve( u, j+istart )
    Dc = Dcurve( u, j+istart )
    SumCurve = bcoeffs[j]*Ac + bcoeffs[j+1]*Bc + bcoeffs[j+2]*Cc + bcoeffs[j+3]*Dc
    return SumCurve

def PrintSumpowerCoeffsLists(bcoeffs):
    print "{ "
    for i in range(nn):
        powerCoeffs_unrounded = SumpowerCoeffs(bcoeffs, i )
        powerCoeffs = [ '%.12f' % elem for elem in powerCoeffs_unrounded ]
        print powerCoeffs[0],powerCoeffs[1],powerCoeffs[2],powerCoeffs[3],""
    powerCoeffs = extendLast(bcoeffs)
    print powerCoeffs[0],powerCoeffs[1],powerCoeffs[2],powerCoeffs[3],"\n};"

def StringSumpowerCoeffsLists(bcoeffs):
    coeffs = ""
    for i in range(nn):
        powerCoeffs_unrounded = SumpowerCoeffs(bcoeffs, i )
        powerCoeffs = [ '%.12f' % elem for elem in powerCoeffs_unrounded ]
        coeffs += str(powerCoeffs[0]) +" "+ str(powerCoeffs[1]) +" "+ str(powerCoeffs[2]) +" "+ str(powerCoeffs[3])+"\n"
    powerCoeffs = extendLast(bcoeffs)
    coeffs += str(powerCoeffs[0]) +" "+ str(powerCoeffs[1]) +" "+ str(powerCoeffs[2]) +" "+ str(powerCoeffs[3])
    return coeffs

# This uses the power coefficients in eavh interval j, to plot the shape of the acceptance.
def calcAcceptance(bcoeffs):
    npt=10
    x = []
    y= []
    for j in range(nn):
        ustart = float(knots[j+istart])
        uend = float(knots[j+istart+1])
        du = (uend-ustart)/npt
        pwr = SumpowerCoeffs(bcoeffs,j)
        for pt in range(npt):
            u = ustart+pt*du
            x.append(u)
            y.append( pwr[0] + pwr[1]*u + pwr[2]*(u**2) + pwr[3]*(u**3)  )

    # artificially extend
    ustart = float(knots[nn+istart])
    uend = 15.
    du = (uend-ustart)/npt
    pwr = extendLast(bcoeffs)
    for pt in range(npt):
        u = ustart+pt*du
        x.append(u)
        y.append( pwr[0] + pwr[1]*u + pwr[2]*(u**2) + pwr[3]*(u**3)  )

    return x,y

# This does the same as the above but uses the formula for S(u) directly to plot the shape of the curve.
def calcCurve(bcoeffs):
    npt=10
    x = []
    y= []
    for j in range(nn):
        ustart = float(knots[j+istart])
        uend = float(knots[j+istart+1])
        du = (uend-ustart)/npt
        for pt in range(npt):
            u = ustart+pt*du
            c = SumCurve(bcoeffs, u, j)
            x.append(u)
            y.append( c )
    return x,y


def extendLast(bcoeffs):
    pwr = SumpowerCoeffs(bcoeffs,nn-1)
    dpwr = [ pwr[1], 2.*pwr[2], 3.*pwr[3], 0. ]
    uend = knots[nn+istart]
    v0 = pwr[0] + pwr[1]*uend + pwr[2]*(uend**2) + pwr[3]*(uend**3)
    m0 = dpwr[0] + dpwr[1]*uend +dpwr[2]*(uend**2)

    coeffs = [ v0-m0*uend, m0, 0., 0. ]
    return coeffs

# This does the same as the above but uses the formula for S(u) directly to plot the shape of the curve.
def main(bcoeffs, inp, trigger):
    #powerCoeffs = SumpowerCoeffs( 0 )
    #print "Petes  seg1: ", powerCoeffs
    f = open("TimeAccCoeff_{}_data_015_{}.txt".format(inp, trigger),"w")
    print "Creating TimeAccCoeff_{}_data_015_{}.txt".format(inp, trigger)
    f.write("0.3\n0.58\n0.91\n1.35\n1.96\n3.01\n7.\n15.\n"+ StringSumpowerCoeffsLists(bcoeffs))
    f.close()

    #PrintSumpowerCoeffsLists(bcoeffs)

    x1,y1 = calcAcceptance(bcoeffs)

    x2,y2 = calcCurve(bcoeffs)
    xx2 = np.array(x2)
    yy2 = np.array(y2) +0.0 #offset just so that it shows as a different line
    pl.plot(x1,y1)
    pl.plot(xx2,yy2)
    pl.axis([0., 15., 0.,2.5])
    #pl.show()

    data = np.loadtxt("hist200.txt")
    #print data
    x=[]
    y=[]
    for item in data :
        x.append((item[0]+item[1])/2.)
        y.append(item[2]*14.-0.05)

    pl.plot(x,y)
    #pl.ylim((0.0,1.5))
    pl.show()

## For what datasets get acceptances
inputs = [ "2015","2016"]

for inp in inputs:
    ntuple = data_loaded[inp]["Ntuple"][0]["Name"].replace("_","\_")
    knots = [float(en["Value"]) for en in data_loaded[inp]["TimeAccParameter"]["KnotParameter"]]
    coef1 = [float(en["Value"]) for en in data_loaded[inp]["TimeAccParameter"]["SplineBiased"]] 
    coef2 = [float(en["Value"]) for en in data_loaded[inp]["TimeAccParameter"]["SplineUnbiased"]] 
    knots.insert(0, knots[0])
    knots.insert(0, knots[0])
    knots.insert(0, knots[0])
    knots.extend([knots[-1],knots[-1],knots[-1]])
    nn=len(coef1)-3  # number of intervals n = len(bcoeffs)-3
    istart=3   # offset into knots list to n=0 knot
    main(coef1, inp, "biased")
    main(coef2, inp, "unbiased")

