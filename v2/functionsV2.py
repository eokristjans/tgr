# -*- coding: utf-8 -*-
"""
Functions used in v2
@author: Erling Oskar
"""
import numpy as np


""" SA 1 """
# tol gives correct to 8 decimal places+
# integral from t1=0 to t2=T
def ComputeArcLengthTPR(f, T, t1=0, tol = 0.5e-8):
    return _AdaptiveIntegrationTPR(f, t1, T, tol)

# Just a function that calls another function
def _AdaptiveIntegrationTPR(f, a, b, tol):
    return _AdQuadTrapezoid(f, f(a), f(b), a, b, tol)

# Computes the AdaptiveQuadrature
def _AdQuadTrapezoid(f, fa, fb, a, b, tol):
    c = (a+b)/2
    h = b-a
    fc = f(c)
    S_ab = (h/2)*(fa+fb)
    S_ac = (h/4)*(fa+fc)
    S_cb = (h/4)*(fc+fb)
    if np.abs(S_ab - S_cb - S_ac) < 3*tol: #  and h < 1e-2
        return (S_ac + S_cb)
    else:
        return ( _AdQuadTrapezoid(f, fa, fc, a, c, tol/2) +
                _AdQuadTrapezoid(f, fc, fb, c, b, tol/2) )
            
        
    
# Possible useful
def ThreePointCenteredDifferenceFormula(f, a=0, h=1e-10):
    return (f(a+h) - f(a-h))/(2*h)



""" SA 2 """
# Usage:    BisectionMethod(f,a,b)
# Pre:      f is a continuous fucntion on the interval [a,b]
#           a < b are numbers such that f(a)*f(b) < 0
#           tol is a number
# Post:     c is the root of f, i.e f(c) = 0
def BisectionMethod(f,a,b,tol=0.5e-6):
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



""" SA 3 """
def ComputeArcLengthSR(f, T, t1=0, tol = 0.5e-8):
    return _AdaptiveIntegrationSR(f, t1, T, tol)

def _AdaptiveIntegrationSR(f, a, b, tol):
    fa = f(a)
    fb = f(b)
    c, fc, sab = _SimpsonsRule(f, a, fa, b, fb)
    return _AdQuadSimpson(f, a, fa, b, fb, tol, sab, c, fc)

def _SimpsonsRule(f, a, fa, b, fb):
    h = (b-a)/6
    c = (a+b)/2
    fc = f(c)
    return (c, fc, h * (fa + 4 * fc + fb))

def _AdQuadSimpson(f, a, fa, b, fb, tol, sab, c, fc):
    h = (b-a)/6
    lc, flc, sac  = _SimpsonsRule(f, a, fa, c, fc)
    rc, frc, scb = _SimpsonsRule(f, c, fc, b, fb)
    if abs(sab - sac - scb) < 15*tol: #  and h < 1e-1
        return sac + scb
    return (_AdQuadSimpson(f, a, fa, c, fc, tol/2, sac , lc, flc) +
           _AdQuadSimpson(f, c, fc, b, fb, tol/2, scb, rc, frc))
    


""" SA 4 """
#
def NewtonsMethod(f, df, xold, tol=0.5e-4):
    xnew = xold - f(xold)/df(xold)
    while abs(xnew-xold) > tol:
        xold = xnew
        xnew = xold - f(xold)/df(xold)
    return xnew



""" Til að prófa Newton. Virkar MJÖG hratt.
fx = lambda x: x**3 + x - 1
dfdx = lambda x: 3*x**2 + 1
ans = NewtonsMethod(fx, dfdx, -0.7)
"""






