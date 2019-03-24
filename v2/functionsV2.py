# -*- coding: utf-8 -*-
"""
Functions used in v2
@author: Erling Oskar
"""
# To install autograd run pip install autograd
# Note: pylint incorrectly gives errors on 
# the autograd.numpy module and grad function
import autograd.numpy as np
from autograd import *

""" SA 1 """
# Returns the lamda function sqrt(dxdt(t)^2 + dydt(t)^2)
# given x(t) and y(t)
def sqrtFunSquared(x, y):
    dxdt = grad(x)
    dydt = grad(y)
    return lambda t: np.sqrt(dxdt(t)^2 + dydt(t)^2)

# Computes the arc length for a given function (f)
# and tolerance (tol), from t = 0 to t = T
# with a given method of integration. 
def compArcLength(f, T, intMethod, tol = 0.5e-8):
    return intMethod(f, 0, T, tol)

# Computes the definite integral of a given function (f)
# from a to b using the method of Adaptive Quadrature
def adQuad(f, a, b, tol):
    c = (a+b)/2
    h = b-a
    S_ab = (h/2)*(f(a)+f(b))
    S_ac = (h/4)*(f(a)+f(c))
    S_cb = (h/4)*(f(c)+f(b))
    if np.abs(S_ab - S_cb - S_ac) < 3*tol: #  and h < 1e-2
        return (S_ac + S_cb)
    else:
        return (adQuad(f, a, c, tol/2) +
                adQuad(f, c, b, tol/2))

""" SA 2 """
# Usage:    BisectionMethod(f,a,b)
# Pre:      f is a continuous fucntion on the interval [a,b]
#           a < b are numbers such that f(a)*f(b) < 0
#           tol is a number
# Post:     c is the root of f, i.e f(c) = 0
def bisectionMethod(f,a,b,tol=0.5e-6):
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

def tStarOfS(f, s, intMethod, rootMethod):
    a = compArcLength(f, 1, intMethod)
    g = lambda b: (s * a - compArcLength(f, b, intMethod))
    return rootMethod(g, 0, 1)

""" SA 3 """
def computeArcLengthSR(f, T, t1=0, tol = 0.5e-8):
    return adQuadSimpson(f, t1, T, tol)

def simpsonsRule(f, a, b):
    h = (b-a)/6
    c = (a+b)/2
    return (c, h * (f(a) + 4 * f(c) + f(b)))

def _adaptiveIntegrationSR(f, a, b, tol):
    
    return adQuadSimpson(f, a, fa, b, fb, tol, sab, c, fc)

def adQuadSimpson(f, a, b, tol, c):
    c, sab = simpsonsRule(f, a, b)
    h = (b-a)/6
    lc, flc, sac  = simpsonsRule(f, a, c)
    rc, frc, scb = simpsonsRule(f, c, b)
    if abs(sab - sac - scb) < 15*tol: #  and h < 1e-1
        return sac + scb
    return (adQuadSimpson(f, a, c, tol, sac, lc, flc) +
           adQuadSimpson(f, c, b, tol, scb, rc, frc))

    
""" SA 4 """
def newtonsMethod(f, xold, tol=0.5e-4):
    dfdt = grad(f)
    xnew = xold - f(xold)/dfdt(xold)
    while abs(xnew-xold) > tol:
        xold = xnew
        xnew = xold - f(xold)/dfdt(xold)
    return xnew