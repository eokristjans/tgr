# -*- coding: utf-8 -*-
"""
Functions used in v2
@author: Erling Oskar
"""

""" Global Imports """
# Plots and animation
import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation, rc
from celluloid import Camera
import matplotlib.pyplot as plt

from time import perf_counter # For timing
import numpy as np
from random import uniform

from sympy.core import sympify
from sympy.utilities import lambdify
from sympy import symbols, diff


    
""" Plots many sections of a curve on one plot """
def plotParameterizedCurves(x, y, ts, n, header, save=None):
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
    if save is not None:
        fig.savefig('img/'+save+'.PNG')
        plt.close()

  

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



""" SA 5 """
# Animates a particle traveling along the path represented by
# function x(t) and y(t)
def animateBezierCurve(x,y,ts,header):
    rc('animation', html='jshtml')
    fig = plt.figure()
    axis = fig.gca()
    axis.set(xlabel='x(t)', ylabel='y(t)',title=header,)
    axis.set_xticks(np.arange(-1.5,2.0,0.5))
    axis.set_yticks(np.arange(-0.5,2.5,0.5))
    plt.close()
    cam = Camera(fig)
    zeroToOne = np.arange(0, 1.01, 0.01)
    for t in ts:
        axis.plot(x(zeroToOne), y(zeroToOne),color='red')
        axis.plot(x(t),y(t), 'bo', color='blue')
        axis.axhline(y=0, color='k')
        axis.axvline(x=0, color='k')
        cam.snap()
    anim = cam.animate(blit=False, interval=100)
    return anim




""" SA 6 """
# Points: List of points (x,y)
# Initial point is: x0 = Points[0][0] and y0 = Points[0][1]
# End point is: x3 = Points[3][0] and y3 = [3][1]
# Control points in between
# Also returns string representation of them
def bezierCurve(Points):
    bx = round(3*(Points[1][0] - Points[0][0]),2)
    cx = round(3*(Points[2][0] - Points[1][0]) - bx,2)
    dx = round(Points[3][0] - Points[0][0] - bx - cx,2)
    by = round(3*(Points[1][1] - Points[0][1]),2)
    cy = round(3*(Points[2][1] - Points[1][1]) - by,2)
    dy = round(Points[3][1] - Points[0][1] - by - cy,2)
    x = lambda t: Points[0][0] + bx*t + cx*t**2 + dx*t**3
    y = lambda t: Points[0][1] + by*t + cy*t**2 + dy*t**3
    xExpr = sympify(str(Points[0][0]) + "+" + str(bx) + "*t" + 
                    "+"+ str(cx) + "*t**2" + "+" + str(dx) + "*t**3")
    yExpr = sympify(str(Points[0][1]) + "+" + str(by) + "*t" + 
                    "+"+ str(cy) + "*t**2" + "+" + str(dy) + "*t**3")
    tSymbol = symbols('t')
    dxdt = lambdify(tSymbol, diff(xExpr, tSymbol), "math")
    dydt = lambdify(tSymbol, diff(yExpr, tSymbol), "math")
    return (x, y), (dxdt, dydt)

# Generates random bezier curves
def bezierCurves(n_BezierCurves):
    BezierCurves = []
    DerivBezierCurves = []
    for i in range(n_BezierCurves):
        points = []
        for j in range(4):
            points.append([round(uniform(-2.0,1.5),2), round(uniform(-1.0,2.0))])
        bc, dbc = bezierCurve(points)
        BezierCurves.append(bc)
        DerivBezierCurves.append(dbc)
        plotParameterizedCurves(BezierCurves[i][0],BezierCurves[i][1],[0,1],1,
                        "Random Bezier Curve "+str(i), "RandomBezierCurve"+str(i))
    return BezierCurves, DerivBezierCurves





""" SA 7 """