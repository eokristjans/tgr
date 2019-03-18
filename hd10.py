import math

""" Dæmi 3 
Adaptive Quadrature as given in Sauer uses the Trapezoid rule to define 
S_[a,b] = 0.5*(b-a)(f(a)+f(b))
Here we will use Simpson's rule
Look better into pages 269-270 in the textbook

I guess we'll basically need 


"""
def SimpsonAdaptiveQuadrature(f, x0, x2, tol):
    x1 = 0.5*(x0+x2)

""" Dæmi 1 """
def RombergIntegration(f, a, b, n):
    R = [[]]
    h = (b-a)/2
    r = h*(f(a)+f(b))
    R[0].append(r)
    for j in range(1,n):
        R.append([])
        h = (b-a)*(0.5**j)
        s = 0
        for i in range(1,2**(j-1)+1):
            s += f(a+(2*i-1)*h)
        r = 0.5*R[j-1][0] + h*s
        R[j].append(r)
        for k in range(1,j+1):
            r = (4**k)*R[j][k-1]-R[j-1][k-1]
            r /= ((4**k)-1)
            R[j].append(r)
    return R
    
ff = lambda x: x*math.exp(x)
aa = 0
bb = 1
Ans = RombergIntegration(ff, aa, bb, 3)
for A in Ans:
    print(A)
