# -*- coding: utf-8 -*-
from functionsV2 import *
import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

# Temporary global tol
tol = 0.5e-5
# Temporary global number of decimals
dec = 5

""" Plot of the given parametric path """
# Given parametric path
x = lambda t: 0.5 + 0.3*t + 3.9*t**2 - 4.7*t**3
y = lambda t: 1.5 + 0.3*t + 0.9*t**2 - 2.7*t**3

plotParameterizedCurves(x,y,[0,1],1,"Original Parametrized Curve")


""" Suggested Activity 1 """
# Manually computed derivatives of the given parametric path
dxdt = lambda t: 0.3 + 7.8*t - 14.1*t**2
dydt = lambda t: 0.3 + 1.8*t - 8.1*t**2

# Compute the lambda function that
# is to be integrated
f = sqrtFunSquared(dxdt, dydt)

# Compute the corresponding arc length
# by computing the integral of f from
# 0 to 1 using the method of Adaptive Quadrature
arcLength = compArcLength(f, 0.0, 1.0, adQuad, tol)

print("The arc length is", round(arcLength, dec))

""" Suggested Activity 2 """
# Let's time the function tStarOfS for s = 0.5 
s = 0.5
start = perf_counter()
tStar2 = tStarOfSBisect(f, s, adQuad, tol)
end = perf_counter()
elapsedTime2 = end - start

print('The optimal value of t for s =', round(s, dec), 'is', round(tStar2, dec))
print('and was computed in', round(elapsedTime2, dec), 'seconds')

# Let's verify if tStar really is the root of the function 
print('Which is the root within our tolerance:', (np.abs(compArcLength(f, 0.0, tStar2, adQuad, tol) / 
                           arcLength - s) < 2*tol))

""" Suggested Activity 3 """
# For n = 4
sArray_n4 = [0.0, 0.25, 0.5, 0.75, 1.0]
n = 4
tStarArray = np.zeros(n+1)
start = perf_counter()
for i in range(n):
    tStarArray[i+1] = tStarOfSBisect(f, sArray_n4[i+1], adQuadSimpson, tol)
end = perf_counter()
elapsedTime3n4 = end - start
print("For n = 4 the computation of t*(s) took",
      round(elapsedTime3n4, dec), "seconds using Simpson's Rule and the Bisection method")

for i in range(n): # let's verify that each arclength is about a quarter of the length of the path
    partialArclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    if abs(partialArclength / arcLength - 1/n) > 10*tol:
        print('Arc length from', round(tStarArray[i], dec), end=' ')
        print('to', round(tStarArray[i+1], dec), end=' ') 
        print('is', round(partialArclength, dec))
        print('Proportional arc length:',
              round(np.abs(partialArclength / arcLength), dec))

plotParameterizedCurves(x,y,tStarArray,n,"Parametrized Curve Equipartitioned into 4 subpaths")
    
# For n = 20
sArray_n20 = np.arange(0.00, 1.05, 0.05)
n = 20
tStarArray = np.zeros(n+1)
start = perf_counter()
for i in range(n):
    tStarArray[i+1] = tStarOfSBisect(f, sArray_n20[i+1], adQuadSimpson, tol)
end = perf_counter()
elapsedTime3n20 = end - start
print("For n = 20 the computation of t*(s) program took", 
      round(elapsedTime3n20, dec), "seconds using Simpson's Rule and the Bisection method")

for i in range(n): # let's verify that each arclength is about one-twentieth of the length of the path
    partialArclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    if abs(partialArclength / arcLength - 1/n) > 10*tol:
        print('Arc length from', round(tStarArray[i], dec), end=' ')
        print('to', round(tStarArray[i+1], dec), end=' ') 
        print('is', round(partialArclength, dec))
        print('Proportional arc length:',
              round(np.abs(partialArclength / arcLength), dec))
    

plotParameterizedCurves(x,y,tStarArray,n,"Parametrized Curve Equipartitioned into 20 subpaths")



""" Suggested Activity 4 """
start = perf_counter() # Newt takes 12 times longer with initial guess of 0.0 than with 0.8
tStar4 = tStarOfSNewton(f, s, adQuad, 0.8, tol)
end = perf_counter()
elapsedTime4 = end - start

print('The computed value of t*(s) for s =', s, end=' ')
print('using Newton\'s method is', round(tStar4, dec+3))

print('The computed value of t*(s) for s =', s, end=' ')
print('using the Bisection method is', round(tStar2, dec+3))
print('Which is the same within our tolerance:', abs(tStar4-tStar2)<tol, "\n")

print('Time required to compute t*(s) using the Newton\'s method is:',
            round(elapsedTime4, dec+2), "seconds")
print('Time required to compute t*(s) using the Bisection method is:',
            round(elapsedTime2, dec+2), "seconds")


# For n = 4
n = 4
tStarArray = np.zeros(n+1)
start = perf_counter()
for i in range(n):
    tStarArray[i+1] = tStarOfSNewton(f, sArray_n4[i+1], adQuadSimpson, tStarArray[i], tol)
end = perf_counter()
elapsedTime4n4 = end - start
print("For n = 4, using Newton's method, the computation of t*(s) took",
      round(elapsedTime4n4, dec), "seconds")
print("For n = 4, using the Bisection method, the computation of t*(s) took",
      round(elapsedTime3n4, dec), "seconds")

for i in range(n): # let's verify that each arclength is about a quarter of the length of the path
    partialArclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    if abs(partialArclength / arcLength - 1/n) > 10*tol:
        print('Arc length from', round(tStarArray[i], dec), end=' ')
        print('to', round(tStarArray[i+1], dec), end=' ') 
        print('is', round(partialArclength, dec))
        print('Proportional arc length:',
              round(np.abs(partialArclength / arcLength), dec))

plotParameterizedCurves(x,y,tStarArray,n,"Parametrized Curve Equipartitioned into 4 subpaths")


# For n = 20
n = 20
tStarArray = np.zeros(n+1)
start = perf_counter()
for i in range(n):
    tStarArray[i+1] = tStarOfSNewton(f, sArray_n20[i+1], adQuadSimpson, tStarArray[i], tol)
end = perf_counter()
elapsedTime4n20 = end - start
print("For n = 20, using Newton's method, the computation of t*(s) took",
      round(elapsedTime4n20, dec), "seconds")
print("For n = 20, using the Bisection method, the computation of t*(s) took",
      round(elapsedTime3n20, dec), "seconds")


for i in range(n): # let's verify that each arclength is about one-twentieth of the length of the path
    partialArclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    if abs(partialArclength / arcLength - 1/n) > 10*tol:
        print('Arc length from', round(tStarArray[i], dec), end=' ')
        print('to', round(tStarArray[i+1], dec), end=' ') 
        print('is', round(partialArclength, dec))
        print('Proportional arc length:',
              round(np.abs(partialArclength / arcLength), dec))

plotParameterizedCurves(x,y,tStarArray,n,"Parametrized Curve Equipartitioned into 20 subpaths")


""" Suggested Activity 5 """
"""
animateBezierCurve(x, y, np.arange(0.0,1.05,0.05),
                   "Traveling at speed with intervals $dt=0.05$ from $t=0$ to $t=1$")
animateBezierCurve(x, y, tStarArray,
                   "Traveling at constant speed given by t*(s)")
"""


""" Suggested Activity 6 """
n_BezierCurves = 3
BezierCurves, DerivBezierCurves = bezierCurves(n_BezierCurves)


n = 20 # Equipartition each BezierCurve into 20 subpaths
tStarArrays = [np.zeros(n+1),np.zeros(n+1),np.zeros(n+1)]
for j in range(n_BezierCurves):
    random_f = sqrtFunSquared(DerivBezierCurves[j][0], DerivBezierCurves[j][1])
    random_arcLength = compArcLength(random_f, 0.0, 1.0, adQuadSimpson, tol)
    print('Arclength of random bezier curve number', j+1, 'is', random_arcLength)
    for i in range(n):
        tStarArrays[j][i+1] = tStarOfSNewton(random_f, sArray_n20[i+1], adQuadSimpson, tStarArrays[j][i], tol)
# Verify that each arclength is about one-twentieth of the length of the path
        partialArclength = compArcLength(random_f, tStarArrays[j][i], tStarArrays[j][i+1], adQuadSimpson, tol)
        if abs(partialArclength / random_arcLength - 1/n) > 10*tol:
            print('Arc length from', round(tStarArrays[j][i], dec), end=' ')
            print('to', round(tStarArrays[j][i+1], dec), end=' ') 
            print('is', round(partialArclength, dec))
            print('Proportional arc length:',
                  round(np.abs(partialArclength / random_arcLength), dec))
    plotParameterizedCurves(BezierCurves[j][0],BezierCurves[j][1], tStarArrays[j], n,
                            "Parametrized Curve Equipartitioned into 20 subpaths",
                            "RandomBezierCurveEquipartitioned"+str(j))

"""
animateBezierCurve(BezierCurves[0][0], BezierCurves[0][1], np.arange(0.0,1.05,0.05),
                 "Traveling at speed with intervals $dt=0.05$ from $t=0$ to $t=1$")
animateBezierCurve(BezierCurves[0][0], BezierCurves[0][1], tStarArrays[0],
                       "Traveling at constant speed given by t*(s)")

animateBezierCurve(BezierCurves[1][0], BezierCurves[1][1], np.arange(0.0,1.05,0.05),
                 "Traveling at speed with intervals $dt=0.05$ from $t=0$ to $t=1$")
animateBezierCurve(BezierCurves[1][0], BezierCurves[1][1], tStarArrays[1],
                       "Traveling at constant speed given by t*(s)")

animateBezierCurve(BezierCurves[2][0], BezierCurves[2][1], np.arange(0.0,1.05,0.05),
                 "Traveling at speed with intervals $dt=0.05$ from $t=0$ to $t=1$")
animateBezierCurve(BezierCurves[2][0], BezierCurves[2][1], tStarArrays[2],
                       "Traveling at constant speed given by t*(s)")

"""




""" Suggested Activity 7 """












