# -*- coding: utf-8 -*-
from functionsV2 import *
import numpy as np
from time import perf_counter

# Temporary global tol
tol = 0.5e-5
# Temporary global number of decimals
dec = 5

""" Suggested Activity 1 """
# Given parametric path with
# x(t) = 0.5 + 0.3*t + 3.9*t**2 - 4.7*t**3 and 
# y(t) = 1.5 + 0.3*t + 0.9*t**2 - 2.7*t**3
# We can manually compute the derivatives
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
# Let's time the function tStarOfS
# for s = 0.5 
s = 0.5
start = perf_counter()
tStar2 = tStarOfSBisect(f, s, adQuad, tol)
end = perf_counter()
elapsedTime2 = end - start

print('The optimal value of t for s=', round(s, dec), ' is ', round(tStar2, dec))
print('and was computed in', round(elapsedTime2, dec), "seconds")

# Let's verify if tStar really is the
# root of the function 
print('Verified:', (np.abs(compArcLength(f, 0.0, tStar2, adQuad, tol) / 
                           compArcLength(f, 0.0, 1.0, adQuad, tol) - s) < 0.001))

""" Suggested Activity 3 """
# For n = 4
sArray = [0.0, 0.25, 0.5, 0.75, 1.0]

start = perf_counter()
tStarArray = [tStarOfSBisect(f, s, adQuadSimpson, tol) for s in sArray]
end = perf_counter()
elapsedTime3n4 = end - start
print("For n = 4 the computation of t*(s) took",
      round(elapsedTime3n4, dec), "seconds")

for i in range(len(sArray)-1):
    arclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    print('Arc length from', round(tStarArray[i], dec), end=' ')
    print('to', round(tStarArray[i+1], dec), end=' ') 
    print('is', round(arclength, dec))
    print('Proportional arc length:', end=' ')
    print(round(np.abs(arclength / compArcLength(f, 0.0, 1.0, adQuadSimpson, tol)), dec))
    
# For n = 20
sArray = np.arange(0.00, 1.05, 0.05)
start = perf_counter()
tStarArray = [tStarOfSBisect(f, s, adQuadSimpson, tol) for s in sArray]
end = perf_counter()
elapsedTime3n20 = end - start
print("For n = 20 the computation of t*(s) program took", 
      round(elapsedTime3n20, dec), "seconds")
 
for i in range(len(sArray)-1):
    arclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    print('Arclength from', round(tStarArray[i], dec), end='')
    print(' to', round(tStarArray[i+1], dec), end=' ') 
    print('is:', round(arclength, dec))
    print('Proportional Arclength :', round(np.abs(arclength / 
                           compArcLength(f, 0.0, 1.0, adQuadSimpson, tol)), 2))
    
    
""" Suggested Activity 4 """
xGuess = 0.8
start = perf_counter()
tStar4 = tStarOfSNewton(f, s, adQuad, xGuess, tol)
end = perf_counter()
elapsedTime4 = end - start

print('The computed value of t*(s) for s =', s, end=' ')
print('using the Newton\'s method is')
print(round(tStar4, dec+3))

print('The computed value of t*(s) for s =', s, end=' ')
print('using the Bisection method is')
print(round(tStar2, dec+3), "\n")

print('Time required to compute t*(s) using the Newton\'s method is:')
print(round(elapsedTime4, dec+2), "seconds")
print('Time required to compute t*(s) using the Bisection method is:')
print(round(elapsedTime2, dec+2), "seconds")
