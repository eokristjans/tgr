import math

""" Dæmi 3 
Adaptive Quadrature as given in Sauer uses the Trapezoid rule to define 
S_[a,b] = 0.5*(b-a)(f(a)+f(b))
Here we will use Simpson's rule
Look better into pages 269-270 in the textbook
Modify the Matlab code for Adaptive Trapezoid Rule Quadrature to use Simpson’s Rule
instead, applying the criterion (5.42) with the 15 replaced by 10. Approximate the integral in
Example 5.12 within 0.005, and compare with Figure 5.5(b). How many subintervals were
required?

|S[a,b] − (S[a,c] + S[c,b])| < 15 ∗ TOL
"""
# Just a function that calls another function
def _AdaptiveIntegration(f, a, b, tol):
    fa = f(a)
    fb = f(b)
    c, fc, sab = _SimpsonsRule(f, a, fa, b, fb)
    return _AdQuadSimpson(f, a, fa, b, fb, tol, sab, c, fc)

def _SimpsonsRule(f, a, fa, b, fb):
    h = (b-a)/6
    c = (a+b)/2
    fc = f(c)
    return (c, fc, h * (fa + 4 * fc + fb))

def _AdQuadSimpson(f, a, fa, b, fb, tol, sab, c, fc):
    h = (b-a)/6
    lc, flc, sac  = _SimpsonsRule(f, a, fa, c, fc)
    rc, frc, scb = _SimpsonsRule(f, c, fc, b, fb)
    if abs(sab - sac - scb) < 15*tol and h < 1e-1:
        return sac + scb
    return (_AdQuadSimpson(f, a, fa, c, fc, tol/2, sac , lc, flc) +
           _AdQuadSimpson(f, c, fc, b, fb, tol/2, scb, rc, frc))

# Dæmi um keyrslu
g = lambda x: 1 + math.sin(math.exp(3*x))
ans = _AdaptiveIntegration(g,-1,1,0.005)
print(ans)
# Gefur 2.5010131869874983




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
