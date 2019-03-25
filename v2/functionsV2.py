# -*- coding: utf-8 -*-
"""
Functions used in v2
@author: Erling Oskar
"""
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter


    
""" Plots many sections of a curve on one plot """
def plotParameterizedCurves(x, y, ts, n, header):
    fig, axis = plt.subplots()
    for i in range(n):
        t = np.arange(ts[i], ts[i+1], 0.001)
        axis.plot(x(t), y(t))
    axis.set(xlabel='x(t)', ylabel='y(t)',title=header)
    plt.yticks(np.arange(-0.5,2.5,0.5))
    plt.xticks(np.arange(-1.5,2.0,0.5))
    axis.axhline(y=0, color='k')
    axis.axvline(x=0, color='k')
    axis.grid()
  

""" SA 1 """
# Returns the lamda function sqrt(dxdt(t)^2 + dydt(t)^2)
# given dxdt(t) and dydt(t)
def sqrtFunSquared(dxdt, dydt):
    return lambda t: np.sqrt(dxdt(t)**2 + dydt(t)**2)

# Computes the arc length for a given function (f)
# and tolerance (tol), from t = 0 to t = T
# with a given method of integration. 
def compArcLength(f, t, T, intMethod, tol):
    return intMethod(f, t, T, tol)

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
        return (adQuad(f, a, c, tol/2) + adQuad(f, c, b, tol/2))

# Returns a function that approximates of the derivative of f(x)
# at x = a
def threePointCentDiff(f, h=1e-5):
    return lambda a: (f(a+h) - f(a-h))/(2*h)

""" SA 2 """
# Usage:    BisectionMethod(f,a,b)
# Pre:      f is a continuous fucntion on the interval [a,b]
#           a < b are numbers such that f(a)*f(b) < 0
#           tol is a number
# Post:     c is the root of f, i.e f(c) = 0
def bisectionMethod(f,a,b,tol):
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

# intMethod: Integration Method
def tStarOfSBisect(f, s, intMethod, tol):
    a = compArcLength(f, 0.0, 1.0, intMethod, tol)
    g = lambda b: (s * a - compArcLength(f, 0.0, b, intMethod, tol))
    return bisectionMethod(g, 0.0, 1.0, tol)

""" SA 3 """
def simpsonsRule(f, a, b):
    h = (b-a)/6
    c = (a+b)/2
    return (c, h * (f(a) + 4 * f(c) + f(b)))

def adQuadSimpson(f, a, b, tol, c=0):
    c, sab = simpsonsRule(f, a, b)
    lc, sac  = simpsonsRule(f, a, c)
    rc, scb = simpsonsRule(f, c, b)
    if abs(sab - sac - scb) < 15*tol: #  and h < 1e-1
        return sac + scb
    return (adQuadSimpson(f, a, c, tol, lc) +
           adQuadSimpson(f, c, b, tol, rc))
    
""" SA 4 """

def newtonsMethod(f, xold, tol):
#    start = perf_counter()
    dfdt = threePointCentDiff(f)
    xnew = xold - f(xold)/dfdt(xold)
    while abs(xnew-xold) > tol:
        xold = xnew
        xnew = xold - f(xold)/dfdt(xold)
#    end = perf_counter()
#    print('it took', end-start, 'time to compute', xnew)
    return xnew

def tStarOfSNewton(f, s, intMethod, xold, tol):
    a = compArcLength(f, 0.0, 1.0, intMethod, tol)
    g = lambda b: (s * a - compArcLength(f, 0.0, b, intMethod, tol))
    return newtonsMethod(g, xold, tol)




""" SA 6 """
# Points: List of points (x,y)
# Initial point is: x0 = Points[0][0] and y0 = Points[0][1]
# End point is: x3 = Points[3][0] and y3 = [3][1]
# Control points in between
def bezierCurve(Points):
    bx = 3*(Points[1][0] - Points[0][0])
    cx = 3*(Points[2][0] - Points[1][0]) - bx
    dx = Points[3][0] - Points[0][0] - bx - cx
    by = 3*(Points[1][1] - Points[0][1])
    cy = 3*(Points[2][1] - Points[1][1]) - by
    dy = Points[3][1] - Points[0][1] - by - cy
    return (lambda t: Points[0][0] + bx*t + cx*t**2 + dx*t**3, lambda t: Points[0][1] + by*t + cy*t**2 + dy*t**3)



