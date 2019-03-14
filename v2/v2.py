# -*- coding: utf-8 -*-

"""
Föll fyrir Verkefni 2
"""
import numpy as np

# tol gives correct to 8 decimal places+
# integral from t1=0 to t2=T
def ComputeArcLength(f, T, t1=0, tol = 0.5e-8):
    return AdaptiveIntegration(f, t1, T, tol)

# Just a function that calls another function
def AdaptiveIntegration(f, a, b, tol):
    return AdaptiveQuadratureS(f, f(a), f(b), a, b, tol)

# Computes the AdaptiveQuadrature
def AdaptiveQuadratureS(f, fa, fb, a, b, tol):
    c = (a+b)/2
    h = b-a
    fc = f(c)
    S_ab = (h/2)*(fa+fb)
    S_ac = (h/4)*(fa+fc)
    S_cb = (h/4)*(fc+fb)
    if np.abs(S_ab - S_cb - S_ac) < 3*tol and h < 1e-2:
        return (S_ac + S_cb)
    else:
        return ( AdaptiveQuadratureS(f, fa, fc, a, c, tol/2) +
                AdaptiveQuadratureS(f, fc, fb, c, b, tol/2) )
            

def SqrtOfFunctionsSquared(x,y):
    return lambda t: np.sqrt(x(t)**2 + y(t)**2)
    

# Possible useful
def ThreePointCenteredDifferenceFormula(f, a=0, h=1e-10):
    return (f(a+h) - f(a-h))/(2*h)

# example functions and usage of computeArcLength
dxdt = lambda t: 0.3 + 7.8*t - 14.1*t**2
dydt = lambda t: 0.3 + 1.8*t - 8.1*t**2
fsa1 = SqrtOfFunctionsSquared(dxdt, dydt)

SA1 = ComputeArcLength(fsa1, 1)
SA1





""" Suggested Activity 2 """

# Notkun:   BisectionMethod(f,a,b)
# Fyrir:    f is a continuous fucntion on the interval [a,b]
#           a < b are numbers such that f(a)*f(b) < 0
#           tol is a number
# Eftir:    c is the root of f, i.e f(c) = 0
def BisectionMethod(f,a,b,tol=0.5e-4):
    fa=f(a)
    fb=f(b)
    while (b-a)/2 > tol:
        c=(a+b)/2
        fc=f(c)
        if fa*fc <= 0:
            b=c
            fb=fc
        else:
            a=c
            fa=fc
    return c

def tStarOfS(f, s):
    arcLen0to1 = ComputeArcLength(f, 1)
    g = lambda b: (s * arcLen0to1 - ComputeArcLength(f, b))
    return BisectionMethod(g, 0, 1)


s = 0.5
ts = tStarOfS(fsa1, s)
print('The optimal value of t for ' + str(s) + ' is ' + str(ts))

# We can verify that with:
print('Verified:', (np.abs(ComputeArcLength(fsa1, ts) / ComputeArcLength(fsa1, 1) - s) < 0.001) )

