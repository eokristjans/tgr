# -*- coding: utf-8 -*-
from functionsV2 import *
import autograd.numpy as np
import time, timeit

""" Suggested Activity 1 """
# Given parametric path with
# x - coordinates x(t) and 
# y - coordinates y(t)
x = lambda t: 0.5 + 0.3*t + 3.9*t**2 - 4.7*t**3
y = lambda t: 1.5 + 0.3*t + 0.9*t**2 - 2.7*t**3

# Compute the lambda function that
# is to be integrgated
f = sqrtFunSquared(x, y)
print(f(1.0))
'''
# Compute the corresponding arc length
# by computing the integral of f from
# 0 to 1 using the method of Adaptive Quadrature
arcLength = compArcLength(f, 1.0, adQuad)

print("The arc length is", arcLength)

""" Suggested Activity 2 """
# Let's time the function tStarOfS
# for s = 0.5 
s = 0.5
start = time.perf_counter()
tStar = tStarOfS(f, s, adQuad, bisectionMethod)
end = time.perf_counter()
time = end - start

print('The optimal value of t for ' + str(s) + ' is ' + str(tStar))
print('and was computed in', time, " seconds")

# Let's verify if tStar really is the
# root of the function 
print('Verified:', (np.abs(compArcLength(f, tStar, adQuad) / 
                           compArcLength(f, 1, adQuad) - s) < 0.001))

""" Suggested Activity 3 """
# For n = 4
sArray = [0.0, 0.25, 0.5, 0.75, 1.0]

# Time the
start = time.perf_counter()
tStarArray = [tStarOfS(f, s, adQuadSimpson, bisectionMethod) for s in sArray]
end = time.perf_counter()
time = end - start 

for i in range(len(sArray)):
    arclength = compArcLength(f, tStarArray[i+1], tStarArray[i])
    print('Arclength from', round(tStarArray[i]), end='')
    print(' to ', round(tStarArray[i+1],2), end='') 
    print('is:', round(arclength,2), end='')
    print('Proportional arc length :', end='')
    print(round(np.abs(arclength / compArcLength(f, 1, adQuadSimpson)), 2))

# For n = 20
sArray = np.arange(0.00, 1.05, 0.05)
time_SA3_n20 = time.perf_counter()
ts3n20 = [tStarOfS_SR(f_SA1, s) for s in s3n20]
time_SA3_n20 = time.perf_counter() - time_SA3_n20

print("SA#, n = 20")
for i in range(20):
    arclength = ComputeArcLengthSR(f_SA1, ts3n20[i+1], ts3n20[i])
    print('Arclength from', round(ts3n20[i], 2), ' to ', round(ts3n20[i+1], 2), 'is:', round(arclength,2))
    print('Proportional Arclength :', round(np.abs(arclength / 
                           ComputeArcLengthSR(f_SA1, 1)), 2))
    
""" Suggested Activity 4 """
def tStarOfS_SR_Newt(f, s):
    a = ComputeArcLengthSR(f, 1)
    g = lambda b: (s * a - ComputeArcLengthSR(f, b))
    return NewtonsMethod(g, 0, 1)

time0_SA4 = time.perf_counter()
ts2_SA4 = tStarOfS_SR_Newt(f_SA1, s2)
time0_SA4 = time.perf_counter() - time0_SA4

print(ts2_SA4)
print(ts2)
print('The optimal value of t for ',s2, 
      ' is equal to the one from Suggest Activity 4:', (abs(ts2_SA4-ts2)<1e-5))
print('Time required to compute t*(s) using the Trapezoid Rule and Newton\'s method',
      'is less than with the Bisection Method:', (time0_SA4 < time_SA2) )



time_SA2
time0_SA4

ts2_SA4
ts2
'''