# -*- coding: utf-8 -*-
from functionsV2 import *
import numpy as np
from time import perf_counter

# Temporary global tol
tol = 0.5e-5

""" Suggested Activity 1 """
# Given parametric path with
# x - coordinates x(t) and 
# y - coordinates y(t)
x = lambda t: 0.5 + 0.3*t + 3.9*t**2 - 4.7*t**3
y = lambda t: 1.5 + 0.3*t + 0.9*t**2 - 2.7*t**3

# Compute the lambda function that
# is to be integrgated
f = sqrtFunSquared(x, y, threePointCentDiff)

# Compute the corresponding arc length
# by computing the integral of f from
# 0 to 1 using the method of Adaptive Quadrature
arcLength = compArcLength(f, 0.0, 1.0, adQuad, tol)

print("The arc length is", round(arcLength, 5))

""" Suggested Activity 2 """
# Let's time the function tStarOfS
# for s = 0.5 
s = 0.5
start = perf_counter()
tStar2 = tStarOfS(f, s, adQuad, bisectionMethod, tol)
end = perf_counter()
elapsedTime2 = end - start

print('The optimal value of t for ' + str(s) + ' is ' + str(tStar2))
print('and was computed in', round(elapsedTime2, 5), "seconds")

# Let's verify if tStar really is the
# root of the function 
print('Verified:', (np.abs(compArcLength(f, 0.0, tStar2, adQuad, tol) / 
                           compArcLength(f, 0.0, 1.0, adQuad, tol) - s) < 0.001))

""" Suggested Activity 3 """
# For n = 4
sArray = [0.0, 0.25, 0.5, 0.75, 1.0]

start = perf_counter()
tStarArray = [tStarOfS(f, s, adQuadSimpson, bisectionMethod, tol) for s in sArray]
end = perf_counter()
elapsedTime3n4 = end - start
print("Í SA3 fyrir n = 4 þá tók keyrslan ", elapsedTime3n4, " sekúndur")

for i in range(len(sArray)-1):
    arclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    print('Arclength from', round(tStarArray[i]), end='')
    print(' to ', round(tStarArray[i+1],2), end='') 
    print('is ', round(arclength,2))
    print('Proportional arc length :', end='')
    print(round(np.abs(arclength / compArcLength(f, 0.0, 1.0, adQuadSimpson, tol)), 2))

# For n = 20
sArray = np.arange(0.00, 1.05, 0.05)
start = perf_counter()
tStarArray = [tStarOfS(f, s, adQuadSimpson, bisectionMethod, tol) for s in sArray]
end = perf_counter()
elapsedTime3n20 = end - start
print("Í SA3 fyrir n = 20 þá tók keyrslan ", elapsedTime3n20, " sekúndur")
 
for i in range(len(sArray)-1):
    arclength = compArcLength(f, tStarArray[i], tStarArray[i+1], adQuadSimpson, tol)
    print('Arclength from', round(tStarArray[i], 2), ' to ', round(tStarArray[i+1], 2), 'is:', round(arclength,2))
    print('Proportional Arclength :', round(np.abs(arclength / 
                           compArcLength(f, 0.0, 1.0, adQuadSimpson, tol)), 2))
    
""" Suggested Activity 4 """
'''
start = perf_counter()
tStar4 = tStarOfS(f, s, adQuad, newtonsMethod, tol)
end = perf_counter()
elapsedTime4 = end - start

print("Í SA4 þá tók keyrslan ", elapsedTime4, " sekúndur")

print('The optimal value of t for ', s, 
      ' is equal to the one from Suggest Activity 4:', (abs(tStar4-tStar2)<1e-5))
print('Time required to compute t*(s) using the Trapezoid Rule and Newton\'s method',
      'is less than with the Bisection Method:', (elapsedTime4 < elapsedTime2) )
'''
