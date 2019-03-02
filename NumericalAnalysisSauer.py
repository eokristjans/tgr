# -*- coding: utf-8 -*-
"""
Projects, functions and notes 
from Numerical Analysis 2nd ed, Sauer, Chapter 4 onwards


@author: Erling Oskar
"""

import numpy as np 
import math
import numpy.linalg as lin
import matplotlib.pyplot as plt
# from fractions import Fraction

'''''''''''''''''''''''''''''''''  Chapter 4  '''''''''''''''''''''''''''''''''
''' Implementation of Gauss-Newton Method which finds the point for which 
the sum of the squared distances to circles is minimized.                   '''
def gaussNewton(r,Dr):
    x = np.zeros(2)
    for k in range(10):
        A = makeA(x,Dr)
        rx = makeR(x,r)
        At = np.transpose(A)
        B = At @ A
        b = -At @ rx
        v = GaussSeidelItr(B,b,10) # Use any method to solve for b
        x += v
    return x

''' 
Params:     r is a n-by-1 array of functions that take 2 arguments each. 
            x=[x0,x1] is a tuple
            applies each function of r to x0 and x1
Returns:    A= [r[0](x0,x1), r[1](x0,x1), ..., r[n](x0,x1)]
'''
def makeR(x, r):
    A = []
    for i in range(len(r)):
        ri=r[i](x[0],x[1])
        A.append([ri])
    A = np.array(A,dtype=float)
    return A
        
''' 
Params:     Dr is a n-by-2 array of functions that take 2 arguments each. 
            x=[x0,x1] is a tuple
            applies each function of Dr to x0 and x1
Returns: A= [[Dr[0][0](x0,x1), Dr[0][1](x0,x1), ...],
             [Dr[1][0](x0,x1), Dr[1][1](x0,x1), ...],
             ... 
             [Dr[n][0](x0,x1), Dr[n][1](x0,x1), ...]]
'''
def makeA(x, Dr):
    A = []
    for i in range(len(Dr)):
        drdxi=Dr[i][0](x[0],x[1])
        drdyi=Dr[i][1](x[1],x[1])
        A.append([drdxi,drdyi])
    A = np.array(A,dtype=float)
    return A

''' A function and its partial derivatives for this assignment '''
def ri(xi,yi,Ri):
    return lambda x,y : np.sqrt((x-xi)**2+(y-yi)**2)-Ri

def drdx(xi,yi):
    return lambda x,y : ((x-(xi)) / (np.sqrt((x-(xi))**2+(y-(yi))**2)))

def drdy(xi,yi):
    return lambda x,y : ((y-(yi)) / (np.sqrt((x-(xi))**2+(y-(yi))**2)))

''' (a)-liður '''
X = [0,1,0]
Y = [1,1,-1]
CircleRadii = [1, 1, 1]

R = []
for i in range(len(CircleRadii)):
    R.append(ri(X[i],Y[i],CircleRadii[i]))
    
DR = []
for i in range(len(CircleRadii)):
    DR.append([drdx(X[i],Y[i]), drdy(X[i],Y[i])])

a = gaussNewton(R,DR)
print(a)


''' (b)-liður '''
X = [-1,1,1]
Y = [0,1,-1]
CircleRadii = [1, 1, 1]

R = []
for i in range(len(CircleRadii)):
    R.append(ri(X[i],Y[i],CircleRadii[i]))
    
DR = []
for i in range(len(CircleRadii)):
    DR.append([drdx(X[i],Y[i]), drdy(X[i],Y[i])])

b = gaussNewton(R,DR)
print(b)

''' Output 
[0.41502111 0.0457165 ]     // Correct answer for (a) is (0.4106, 0.0555)
[2.75549395e-01 9.94749856e-18]     // Corr ans for (b) is (0.2755, 0.0)

Being off by 0.01 in (a) may be caused by high condition number due to
multiply A and At together, or from imperfection in GaussSeidelItr
'''






''' Notadi thetta agaeta fall af heimadaemum 4 til adstodar og upprifjunar

# Notkun:   x = GaussSeidelItr(AA,bb,s)
# Fyrir:    AA er nxn fylki, b er vigur að lengd n, s er jákvæð heiltala
# Eftir:    x er s ítranir af Gauss-Seidel með upphafsvigur np.zeros(n)
'''
def GaussSeidelItr(AA,bb,s):
    A=np.array(AA,float)
    b=np.array(bb,float)
    (n,m)=np.shape(A)
    nb=len(b)
    if n != m or m != nb:
        print('wrong dimensions')
        return
    x_old=np.zeros(n)
    x_new=np.zeros(n)
    for j in range(s):
        for i in range(n):
            x_new[i]=b[i]
            for k in range(n):
                if k != i:
                    x_new[i] -= A[i,k]*x_new[k]     # Breytti hér í x_new
            x_new[i] /= A[i,i]
        for i in range(n):
            x_old[i]=x_new[i]
    return x_new






'''     Implementation of QR-factorization using Householder reflectors.    '''
'''     Only works for nxn array. Requires transposation adjustments.       '''
def householder(AA):
    A = np.array(AA, dtype=float)
    R = np.transpose(A)
    n = len(A)
    H = []
    for i in range(n-1):
        R = np.transpose(R)
        x = R[i][i:n]
        print('x', x)
        v = np.zeros(n-i)
        v[0] = lin.norm(x)
        v -= x
        vInner = np.inner(v,v)
        if vInner < 1.0-10 :
            return 'algorithm cancelled to prevent division by zero'
        print('v', v)
        v = np.concatenate((np.zeros(i),v))
        P = np.outer(v,v) / vInner
        print('P')
        print(P)
        P = np.identity(n) - 2*P
        P = changeP(P,i)
        H.append(P)
        print('H[i]')
        print(H[i])
        R = H[i] @ np.transpose(R)
        print('R')
        print(R)
    Q = H[n-2]
    for j in range(n-3,-1,-1):
        Q = H[j] @ Q
    return (Q,R)

    
def changeP(PP,i):
    P = np.copy(PP)
    for j in range(i):
        P[j][j]=1
    return P
    

# A = [[3,1],[4,3]]
# A = np.transpose(A)
A = np.reshape([4,0,3,8,2,6,1,-2,7],(3,3))
QR = householder(A)
print('A')
print(A)
print('Q')
print(QR[0])
print('R')
print(QR[1])
print('Q @ R')
print(QR[0] @ QR[1]) 

''' Outcome:
R
[[ 5.00000000e+00  1.00000000e+01  5.00000000e+00]
 [ 0.00000000e+00  2.00000000e+00 -2.00000000e+00]
 [ 2.22044605e-16  4.44089210e-16  5.00000000e+00]]
A
[[ 4  0  3]
 [ 8  2  6]
 [ 1 -2  7]]
Q
[[ 0.8  0.  -0.6]
 [ 0.   1.   0. ]
 [ 0.6  0.   0.8]]
R
[[ 5.00000000e+00  1.00000000e+01  5.00000000e+00]
 [ 0.00000000e+00  2.00000000e+00 -2.00000000e+00]
 [ 2.22044605e-16  4.44089210e-16  5.00000000e+00]]
Q @ R
[[ 4.  8.  1.]
 [ 0.  2. -2.]
 [ 3.  6.  7.]]
'''


def classicGramSchmidt(AA):
    A = np.array(AA, dtype=float)
    n = len(A)
    R = np.zeros([n,n], dtype=float)
    Q = np.zeros([n,n], dtype=float)
    for j in range(n):
        y = A[j]
        for i in range(j):
            R[i][j] = np.dot(np.transpose(Q[i]),A[j])
            y -= R[i][j]*Q[i]
        R[j][j] = lin.norm(y)
        Q[j] = y / R[j][j]
    return (np.transpose(Q),R)

A = [[4,0,3],[8,2,6],[1,-2,7]]
B = classicGramSchmidt(A)
print(np.transpose(np.array(A, dtype=float)))
print(A[0])
print(B[0])
print(B[1])
print((np.matmul(B[0],B[1])))
B[1][2][2] *=(-1)



def modifiedGramSchmidt(AA):
    A = np.array(AA, dtype=float)
    n = len(A)
    R = np.zeros([n,n], dtype=float)
    Q = np.zeros([n,n], dtype=float)
    for j in range(n):
        y = A[j]
        for i in range(j):
            R[i][j] = np.dot(np.transpose(Q[i]),y)
            y -= R[i][j]*Q[i]
        R[j][j] = lin.norm(y)
        Q[j] = y / R[j][j]
    return (np.transpose(Q),R)

A = [[4,0,3],[8,2,6],[1,-2,7]]
B = modifiedGramSchmidt(A)
print(np.transpose(np.array(A, dtype=float)))
print(B[0])
print(B[1])
print((np.matmul(B[0],B[1])))
