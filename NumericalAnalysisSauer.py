# -*- coding: utf-8 -*-
"""
Projects, functions and notes 
from Numerical Analysis 2nd ed, Sauer, 
Chapter 4 : Least Squares (Aðdferð minnstu kvaðrata), QR-factorization & more
Chapter 5 : Numerical differentiation and integration

@author: Erling Oskar
"""import numpy as np 
import numpy.linalg as lin
import matplotlib.pyplot as plt
# from fractions import Fraction


'''''''''''''''''''''''''''''''''  Chapter 6  '''''''''''''''''''''''''''''''''
# Runge-Kutta 4 Good iterative solver for Ordinary Differential Equations
def RK4(f, t, w, h):
    g = h/2
    s1 = f(t,w)
    s2 = f(t+g, w+g*s1)
    s3 = f(t+g, w+g*s2)
    s4 = f(t+h, w+h*s3)
    return (t+h, w + h*(s1+2*(s2+s3)+s4)/6)
    

f = lambda t,y: t*y + t**3
y0 = 1
h=0.2
t=0

w = y0
for i in range(5):
    t, w = RK4(f, t, w, h)
    print(t, w)
    
""" Output:
0.2 1.0206033466666666
0.4 1.089858524300753
0.6 1.2316459785404794
0.8 1.4913713892289726
1.0 1.9461400240300424
"""


















def ExplicitTrapezoidMethod(f, y0, h, k, t0=0):
    t = t0
    wCurr = y0
    for i in range(k):
        ftw = f(t,wCurr)
        t +=h
        wTilda = wCurr + h*ftw
        wNext = wCurr + 0.5*h*(ftw + f(t,wTilda))
        wCurr = wNext
        print('w['+str(i+1)+'] =', wCurr)
    return wCurr


kk = 4
hh = 0.25

f_1a = lambda t,y: 2*(t+1)*y
y0_1a = 1
ans_1a = ExplicitTrapezoidMethod(f_1a, y0_1a, hh, kk)
""" Output:
w[1] = 1.71875
w[2] = 3.30322265625
w[3] = 7.070960998535156
w[4] = 16.793532371520996
"""

f_1b = lambda t,y: t**3/y**2
y0_1b = 1
ans_1b = ExplicitTrapezoidMethod(f_1b, y0_1b, hh, kk)
""" Output:
w[1] = 1.001953125
w[2] = 1.019342601468375
w[3] = 1.0822649536177262
w[4] = 1.2182419316985125
"""

f_2b = lambda t,x: (2*t-2)*x
y0_2b = 1
ans_2b = ExplicitTrapezoidMethod(f_2b, y0_2b, hh, kk)
""" Output:
w[1] = 0.65625
w[2] = 0.48193359375
w[3] = 0.39910125732421875
w[4] = 0.3741574287414551
"""






'''''''''''''''''''''''''''''''''  Chapter 5  '''''''''''''''''''''''''''''''''
''' Composite Simpson's Rule
'''
def compositeSimpsonsRule(f,d4f,a,b,m):
    h = (b-a)/(2*m)
    x = np.arange(a,b+h,h)
    c = (b-a)/2 # or some other number between a and b
    Oh4 = (b-a)*h**4*d4f(c)/180
    ySum = 0
    for i in range(1,m+1):
        ySum += 4*f(x[2*i-1])
    for i in range(1,m):
        ySum += 2*f(x[2*i])
    ySum = ySum + f(x[0]) + f(x[2*m])
    return ( ySum*h/3 - Oh4 )


''' Composite Trapezoid Rule 
    Estimates the integral of f from a to b
    f is a mathematical function that takes 1 parameter
    d2f is the second derivative of f
    a and b are the endpoints of the integral
    m is the number of panels that [a,b] will be split into
'''
def compositeTrapezoidRule(f,d2f,a,b,m):
    h = (b-a)/m
    x = np.arange(a,b+h,h)
    c = (b-a)/2 # or some other number between a and b
    Oh2 = (b-a)*h*h*d2f(c)/12
    ySum = 0
    for i in range(1,m):
        ySum += f(x[i])
    ySum = ySum*2 + f(x[0]) + f(x[m])
    return (0.5*h*ySum + Oh2)
    


ff = lambda x: (x**2+4)**(-0.5)
d2ff = lambda x: 3*x**2*(x**2+4)**(-2.5)-(x**2+4)**(-1.5)
d4ff = lambda x: ((2*x**2-4)/((x**2+4)**(5/2)))
aa = 0
bb = 2*3**(0.5)
mm = 16

ff0 = lambda x: 0
FF = lambda x: np.arcsinh(x/2)

cTR = compositeTrapezoidRule(ff,d2ff,aa,bb,mm)
cSR = compositeSimpsonsRule(ff,d4ff,aa,bb,mm)

cTR0 = compositeTrapezoidRule(ff,ff0,aa,bb,mm)
cSR0 = compositeSimpsonsRule(ff,ff0,aa,bb,mm)

rettSvar = FF(bb) - FF(aa)
skekkja = [np.abs(rettSvar-cTR), np.abs(rettSvar-cSR),
           np.abs(rettSvar-cTR0), np.abs(rettSvar-cSR0)]

print('Rétt lausn reiknuð með stofnfalli:', rettSvar)
print()
print('Áætlum c = (b-a)/2')
print('Reiknuð með Samsettu Trapissureglunni:', cTR)
print('Skekkjan er:', skekkja[0])
print('Reiknuð með Samsettu Reglu Simpsons:', cSR)
print('Skekkjan er:', skekkja[1])
print()
print('Án áæltaðar skekkju')
print('Reiknuð með Samsettu Trapissureglunni:', cTR0)
print('Skekkjan er:', skekkja[2])
print('Reiknuð með Samsettu Reglu Simpsons:', cSR0)
print('Skekkjan er:', skekkja[3])
print()
print('Minnsta skekkjan er:', min(skekkja), 'sem er aðferð númer',
      (skekkja.index(min(skekkja))+1), 'af þeim sem við prófuðum' )
print()
print('Aðferð Simpson\'s án áætlaðrar skekkju gefur hér nákvæmasta niðurstöðu')


''' Given Program. What does it do? 
Estimates the difference between a numerical estimation of the derivative of 
sin(x), and the actual value of dsin(x)/dx = cos(x) in point a. 
The computer has problems subtracting floating point numbers of similar value 
from each other, and so it may render a very large number.
The computer also has problems dividing by very small numbers, so it may
render very large numbers even if the actual limit as the denominator 
approaches zero is zero.
'''
def hge1(a):
    p=np.arange(-20,1)
    f=lambda x: np.sin(x)
    df=lambda x: np.cos(x)
    err1,err2,err3=[],[],[]
    for i in range(len(p)):
        h=10.0**p[i]
        err1+=[(f(a+h) - f(a-h))/(2*h)-df(a)]
        err2+=[(f(a+h) - f(a))/h-df(a)]
        h=2*h
        err3+=[(f(a-h)-8*f(a-h/2)+8*f(a+h/2) - f(a+h))/(6*h)-df(a)]
    plt.plot(p,np.log10(np.abs(err1)),'b')
    plt.plot(p,np.log10(np.abs(err2)),'r')
    plt.plot(p,np.log10(np.abs(err3)),'g')
    plt.show()


for i in range(-10,10):
    print(i)
    hge1(i)

hge1(-100)


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
