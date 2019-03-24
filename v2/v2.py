# -*- coding: utf-8 -*-
from functionsV2 import *
import numpy as np
import time, timeit

""" Suggested Activity 1 """
# Given parametric path with
# x - coordinates x(t) and 
# y - coordinates y(t)
x = lambda t: 0.5 + 0.3*t + 3.9*t**2 - 4.7*t**3
y = lambda t: 1.5 + 0.3*t + 0.9*t**2 - 2.7*t**3

# Compute the lambda function that
# is to be integrgated
f1 = sqrtFunSquared(x, y)

# Compute the corresponding arc length
# by computing the integral of f1 from
# 0 to 1 using the method of Adaptive Quadrature
arcLength1 = compArcLength(f1, 1, adQuad)

""" Suggested Activity 2 """

s2 = 0.5
start = time.perf_counter()
ts2 = tStarOfS(f1, s2, adQuad, bisectionMethod)
end = time.perf_counter()
time2 = end - start

print('The optimal value of t for ' + str(s2) + ' is ' + str(ts2))

# We can verify that with:
print('Verified:', (np.abs(ComputeArcLengthTPR(f_SA1, ts2) / 
                           ComputeArcLengthTPR(f_SA1, 1) - s2) < 0.001))

""" Suggested Activity 3 """

def tStarOfS_SR(f, s):
    a = ComputeArcLengthSR(f, 1)
    g = lambda b: (s * a - ComputeArcLengthSR(f, b))
    return BisectionMethod(g, 0, 1)

# For n=4
s3n4 = [0.0, 0.25, 0.5, 0.75, 1.0]

time_SA3_n4 = time.perf_counter()
ts3n4 = [tStarOfS_SR(f_SA1, s) for s in s3n4]
time_SA3_n4 = time.perf_counter() - time_SA3_n4

print("SA3, n = 4")
for i in range(4):
    arclength = ComputeArcLengthSR(f_SA1, ts3n4[i+1], ts3n4[i])
    print('Arclength from', ts3n4[i], ' to ', round(ts3n4[i+1],2), 'is:', round(arclength,2))
    print('Proportional Arclength :', round(np.abs(arclength / 
                           ComputeArcLengthSR(f_SA1, 1)), 2))

# For n=20
s3n20 = np.arange(0.00,1.05,0.05)

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
