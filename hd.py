import numpy as np 
import math
import numpy.linalg as lin
import matplotlib.pyplot as plt
# from fractions import Fraction


''''''''''''''''''''''''''''''''' Dæmablað 7 '''''''''''''''''''''''''''''''''
''' 1. Forritið bæði klassíska og betrumbætta QR-þáttun  '''
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
print(B[0])
print(B[1])
print((np.matmul(B[0],B[1])))

''' Úttak
[[ 4.  8.  1.]
 [ 0.  2. -2.]
 [ 3.  6.  7.]]
[[ 0.8  0.  -0.6]
 [ 0.   1.   0. ]
 [ 0.6  0.   0.8]]
[[ 5. 10.  5.]
 [ 0.  2. -2.]
 [ 0.  0.  5.]]
[[ 4.  8.  1.]
 [ 0.  2. -2.]
 [ 3.  6.  7.]]
'''


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

''' Úttak
[[ 4.  8.  1.]
 [ 0.  2. -2.]
 [ 3.  6.  7.]]
[[ 0.8  0.  -0.6]
 [ 0.   1.   0. ]
 [ 0.6  0.   0.8]]
[[ 5. 10.  5.]
 [ 0.  2. -2.]
 [ 0.  0.  5.]]
[[ 4.  8.  1.]
 [ 0.  2. -2.]
 [ 3.  6.  7.]]
'''





''''''''''''''''''''''''''''''''' Dæmablað 6 '''''''''''''''''''''''''''''''''

chebyshev_bruunarpunktar = lambda j,n: np.cos((np.pi*(2*j-1))/(2*n))

d2a = []
for i in range(1,7):
    d2a.append(chebyshev_bruunarpunktar(i,6))


np.cos((np.pi*(2-1))/(12))

4**4/2**7

def helmadf(f,a,b,TOL):
    fa=f(a)
    fb=f(b)
    while (b-a)/2 > TOL:
        c=(a+b)/2
        fc=f(c)
#        print(fc)
        if fa*fc <= 0:
            b=c
            fb=fc
        else:
            a=c
            fa=fc
    return c,fc
    
TOL = 0.5e-6

f = lambda x: 32*x**6 - 48*x**4 + 18*x**2 - 1
t2 = np.arange(-0.8, 1, 0.1)
fig2, ax2 = plt.subplots()
ax2.plot(t2, f(t2))

ax2.grid()
plt.show()

d2a1 = helmadf(f,-9000,-7500,TOL)
d2a2 = helmadf(f,-7000,-6500,TOL)
d2a1 = helmadf(f,-9000,-7500,TOL)

f(0.9659258)


######## Dæmablað 5 #########
''' T.Gr. hd5-1
Ekki alveg að virka rétt held ég... 
Best að leysa þetta bara með aðferð úr stæ3 fyrst, ekki að ég muni þær lol

def myI(n):
    I = []
    for i in range(n):
        I.append([])
        for j in range(n):
            if i==j:
                I[i].append(1.0)
            else:
                I[i].append(0.0)
    return I

Leysa u^3 - v^3 + u = 0, u^2+v^2 = 1 
með x0=[1,1] og A0=B0=I
með Newton og Broyden I og II
Hvað þarf margar ítranir með hverri umferð
til að fá skekkju innan við 10^(-15) ?
'''
import numpy as np
import numpy.linalg as lin

def Broyden2(f,xx0,k):
    x0 = np.array(xx0)
    n = len(x0)
    B = np.eye(n)
    fx0 = np.array(f(x0))
    for i in range(k):
        x=x0-B@fx0
        d = x-x0
        delta = np.array(f(x))-fx0
        Bdelta = B@delta
        dT = np.transpose(d)
        B = B + ((d-Bdelta)*dT*B)/(dT*Bdelta)
        x0 = x
        print(x0)
    return x0

myf = lambda x: [x[0]**3 - x[1]**3 + x[0], x[0]**2 + x[1]**2 -1]
myx0 = [1,1]
xk = Broyden2(myf,myx0,20)


##########  Dæmablað 4  ##########
'''
#####       Dæmi 1
Notið eftirfarandi forrit og breytið því til þess að nota Gauss-Seidel og SOR aðferðina
Notið svo öll þessi forrit til þess að leysa Example 2.24 í bókinni og gangið úr skugga um að þið fáið
sömu niðurstöðu og bókin
'''
# Notkun:   x = matitr(AA,bb,s)
# Fyrir:    AA er nxn fylki, b er vigur að lengd n, s er jákvæð heiltala
# Eftir:    x er s ítranir af Jacobi með upphafsvigur np.zeros(n)
def matitr(AA,bb,s):
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
                    x_new[i] -= A[i,k]*x_old[k]
            x_new[i] /= A[i,i]
        for i in range(n):
            x_old[i]=x_new[i]
    return x_new


 
'''
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



'''
Ekki rétt
'''
def SORItr(AA,bb,s,omega):  # baetti vid omega innatki
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
            x_new[i]=omega*b[i]                     # omega margfeldi her?
            for k in range(n):
                if k != i:
                    x_new[i] -= omega*A[i,k]*x_new[k]           # omega margfeldi her?
#            x_new[i] *= omega                   # Baetti vid omega margfeldinu
            x_new[i] /= A[i,i]
        for i in range(n):
            x_new[i] += ((1-omega)*x_old[i])    # legg vid omega margfeldid
            x_old[i] =  x_new[i]
    return x_new



A=[[3,-1,0,0,0,1/2], [-1,3,-1,0,1/2,0], [0,-1,3,-1,0,0], [0,0,-1,3,-1,0], [ 0,1/2,0,-1,3,-1], [ 1/2,0,0,0,-1,3]]
b=[5/2,3/2,1,1,3/2,5/2]
s=6
omega=1.1


# Prentar úttakið eins og í bókinni
def db4_1(AA,bb,s,omega):
    jL = matitr(A,b,s)
    gsL= GaussSeidelItr(A,b,s)
    sorL = SORItr(A,b,s,omega)
    n = len(jL)
    print("Jacobi Gau-Se Sor")
    for i in range(n):
        print("%.4f %.4f %.4f" % (jL[i],gsL[i],sorL[i]))
    
db4_1(A,b,s,omega)



''' Haldast óbreyttir
print(np.array(A))
print(b)
print(lin.solve(A,b))   # [1. 1. 1. 1. 1. 1.]
'''


''' Lausn úr bók
Jacobi Gauss–Seidel SOR
0.9879 0.9950 0.9989
0.9846 0.9946 0.9993
0.9674 0.9969 1.0004
0.9674 0.9996 1.0009
0.9846 1.0016 1.0009
0.9879 1.0013 1.0004
'''



########## Dæmablað 3 ##########
##### Dæmi 3
'''
 Computer Problems 2.3.1
For the nxn matrix with entries A[i][j] = 5/(i+2(j+1)), set x = [1,...,1] and Ax=b.
Computer the double precision computed solution, x_c.
Find the 
    infinity norm (óendanlega normið) of the forward error and the error magnification factor of the problem Ax=b, 
and compare it with the condition number of A: (a) n=6, (b) n = 10

'''

# Notkun:   (L,U) = LU(A)
# Fyrir:    A er np.array fylki
# Eftir:    L er nedra thrihyrningsfylki
#           U er efra thrihyrningsfylki
#           L*U = A
#   Haettir keyrslu ef forritid rekst a 0-forustustudul
def LU(A):
    n = len(A)
    L = []
    myI(L,n)
    for j in range(n):
        for i in range(j+1,n):
            mDen = A[j][j]
            # Haettir frekar en ad deila med 0
            if mDen == 0.0:
                print("Forrit haetti keyrslu snemma. Getur ekki deilt med 0")
                return L,A
            mNum = A[i][j]
            m = mNum/mDen
            L[i][j]=m
            for k in range(j,n):
                A[i][k] -= m*A[j][k]
    return L,A


def cp2_3_1(A,n):
    for i in range(n):
        A.append([])
        for j in range(n):
            A[i].append(5/(i+2*(j+1)))


A6 = []
cp2_3_1(A6,6)
A6 = np.array(A6)
#print(A6)

print("n = 6 gefur eftirfarnadi") 
print("Ástandstala (condition number) með n=6 er:", lin.cond(A6,math.inf))


x6 = np.ones(6)
#print(x6)

b6 = A6 @ x6
#print(b6)


A6LU = LU(A6)
A6L = np.array(A6LU[0])
A6U = np.array(A6LU[1])
#print(A6L)
#print(A6U)

y6c = lin.solve(A6L,b6)
#print(y6c)

x6c = lin.solve(A6U,y6c)
#print(x6c)

FE6num = lin.norm(x6c - x6, math.inf)
FE6den = lin.norm(x6, math.inf)
FE6 = FE6num/FE6den
print("Forward Error með LU þáttún:", FE6)

BE6num = lin.norm((A6 @ x6c) - b6, math.inf)
BE6den = lin.norm(b6, math.inf)
BE6 = BE6num/BE6den
print("Backward Error með LU þáttun:", BE6)

ErrMagFac6 = FE6 / BE6
print("Error Magnification Factor með LU þáttun:", ErrMagFac6)


A6 = []
cp2_3_1(A6,6)
A6 = np.array(A6)
#print(A6)

x6 = np.ones(6)
#print(x6)

b6 = A6 @ x6
#print(b6)

x6c = lin.solve(A6,b6)
#print(x6c)

FE6num = lin.norm(x6c - x6, math.inf)
FE6den = lin.norm(x6, math.inf)
FE6 = FE6num/FE6den
print("Forward Error án LU þáttun:", FE6)

BE6num = lin.norm((A6 @ x6c) - b6, math.inf)
BE6den = lin.norm(b6, math.inf)
BE6 = BE6num/BE6den
print("Backward Error án LU þáttun:", BE6)

ErrMagFac6 = FE6 / BE6
print("Error Magnification Factor án LU þáttun:", ErrMagFac6)

print()





A10 = []
cp2_3_1(A10,10)
A10 = np.array(A10)
#print(A10)

print("n = 10 gefur eftirfarnadi") 
print("Ástandstala (condition number) með n=10 er:", lin.cond(A10,math.inf))

x10 = np.ones(10)
#print(x10)

b10 = A10 @ x10
#print(b10)


A10LU = LU(A10)
A10L = np.array(A10LU[0])
A10U = np.array(A10LU[1])
#print(A10L)
#print(A10U)

y10c = lin.solve(A10L,b10)
#print(y10c)

x10c = lin.solve(A10U,y10c)
#print(x10c)

FE10num = lin.norm(x10c - x10, math.inf)
FE10den = lin.norm(x10, math.inf)
FE10 = FE10num/FE10den
print("Forward Error með LU þáttún:", FE10)

BE10num = lin.norm((A10 @ x10c) - b10, math.inf)
BE10den = lin.norm(b10, math.inf)
BE10 = BE10num/BE10den
print("Backward Error með LU þáttun:", BE10)

ErrMagFac10 = FE10 / BE10
print("Error Magnification Factor með LU þáttun:", ErrMagFac10)


A10 = []
cp2_3_1(A10,10)
A10 = np.array(A10)
#print(A10)

x10 = np.ones(10)
#print(x10)

b10 = A10 @ x10
#print(b10)

x10c = lin.solve(A10,b10)
#print(x10c)

FE10num = lin.norm(x10c - x10, math.inf)
FE10den = lin.norm(x10, math.inf)
FE10 = FE10num/FE10den
print("Forward Error án LU þáttun:", FE10)

BE10num = lin.norm((A10 @ x10c) - b10, math.inf)
BE10den = lin.norm(b10, math.inf)
BE10 = BE10num/BE10den
print("Backward Error án LU þáttun:", BE10)

ErrMagFac10 = FE10 / BE10
print("Error Magnification Factor án LU þáttun:", ErrMagFac10)

'''
Útkoma:
n = 6 gefur eftirfarnadi
Ástandstala (condition number) með n=6 er: 70342013.93053748
Forward Error með LU þáttún: 1.587638909228417e-10
Backward Error með LU þáttun: 0.6666666666666667
Error Magnification Factor með LU þáttun: 2.381458363842625e-10
Forward Error án LU þáttun: 4.4538306376296077e-10
Backward Error án LU þáttun: 1.450087215836939e-16
Error Magnification Factor án LU þáttun: 3071422.5937500005

n = 10 gefur eftirfarnadi
Ástandstala (condition number) með n=10 er: 131337064464496.02
Forward Error með LU þáttún: 0.00032580806187954003
Backward Error með LU þáttun: 0.6745937551377676
Error Magnification Factor með LU þáttun: 0.00048296928247283064
Forward Error án LU þáttun: 0.00026835943707226306
Backward Error án LU þáttun: 1.2129573866111993e-16
Error Magnification Factor án LU þáttun: 2212439118096.4287
'''



# DrSiggi
def LUsolP(AA,bb):
    A=np.array(AA,float)
    b=np.array(bb,float)
    eps=0.0
    n=len(A)
    p=[]
    for i in range(n):
        p+=[i]
    for j in range(n-1):
        maxval=abs(A[p[j]][j])
        maxindex=j
        for i in range(j+1,n):
            if abs(A[p[i]][j])>maxval:
                maxval=abs(A[p[i]][j])
                maxindex=i
        p[j],p[maxindex]=p[maxindex],p[j]
        #print(A[p[j]][j])
        if abs(A[p[j]][j])<=eps:
            print('(almost) zero pivot ',A[p[j]][j],)
            return
        for i in range(j+1,n):
            m=-A[p[i]][j]/A[p[j]][j]
            A[p[i]][j]=-m
            for k in range(j+1,n):
                A[p[i]][k]+=m*A[p[j]][k]
    y=np.zeros(n)
    for i in range(n):
        for j in range(i):
            b[p[i]]-=A[p[i]][j]*y[j]
        y[i]=b[p[i]]
    x=np.zeros(n)
    for i in range(n-1,-1,-1):
       for j in range(i+1,n):
           y[i]-=A[p[i]][j]*x[j]
       x[i]=y[i]/A[p[i]][i]
    L=np.eye(n,n)
    U=np.zeros((n,n))
    P=np.zeros((n,n))
    for i in range(n):
        for j in  range(i):
            L[i][j]=A[p[i]][j]
        for j in range(i,n):
            U[i][j]=A[p[i]][j]
        P[i][p[i]]=1.0
    return x.reshape((n,1)),p,P,L,U     


##### Dæmi 2
'''
Write a script to take a matrix A as input and output L and U.
No row exchanges. Program shuts down if zero pivot.
'''

# Breytir I=[] i einingafylki af vidd n
# ath ad I skal vera fyrirfram skilgreint I=[] (tomur listi)
def myI(I,n):
    for i in range(n):
        I.append([])
        for j in range(n):
            if i==j:
                I[i].append(1.0)
            else:
                I[i].append(0.0)

# Prentar fylki a lesanlegri mata
def printArray(A):
    n = len(A)
    for i in range(n):
        print(A[i])
    print()
            

# Notkun:   (L,U) = LU(A)
# Fyrir:    A er np.array fylki
# Eftir:    L er nedra thrihyrningsfylki
#           U er efra thrihyrningsfylki
#           L*U = A
#   Haettir keyrslu ef forritid rekst a 0-forustustudul
def LU(A):
    n = len(A)
    L = []
    myI(L,n)
    for j in range(n):
        for i in range(j+1,n):
            mDen = A[j][j]
            # Haettir frekar en ad deila med 0
            if mDen == 0.0:
                print("Forrit haetti keyrslu snemma. Getur ekki deilt med 0")
                return L,A
            mNum = A[i][j]
            m = mNum/mDen
            L[i][j]=m
            for k in range(j,n):
                A[i][k] -= m*A[j][k]
    return L,A
    
# Prufufall 1
A1 = [[4.0,2.0,0.0],[4.0,4.0,2.0],[2.0,2.0,3.0]]
luA1 = LU(A1)
print("L:")
printArray(luA1[0])
print("U:")
printArray(luA1[1])

# Prufufall 2
A2 = [[4.0,2.0,0.0],[4.0,2.0,2.0],[2.0,2.0,3.0]]
luA2 = LU(A2)
print("L:")
printArray(luA2[0])
print("U:")
printArray(luA2[1])

''' UTTAK:
L:
[1.0, 0.0, 0.0]
[1.0, 1.0, 0.0]
[0.5, 0.5, 1.0]

U:
[4.0, 2.0, 0.0]
[0.0, 2.0, 2.0]
[0.0, 0.0, 2.0]

Forrit haetti keyrslu snemma. Getur ekki deilt med 0
L:
[1.0, 0.0, 0.0]
[1.0, 1.0, 0.0]
[0.5, 0.0, 1.0]

U:
[4.0, 2.0, 0.0]
[0.0, 0.0, 2.0]
[0.0, 1.0, 3.0]
'''



########## Dæmablað 2 ##########
##### Dæmi 3

# Notkun:   xn = secant_method(f,x0,x1,n)
# Fyrir:    f er samfelld einundaradgerd
#           skilgreind á  bilinu [a,b].
#           x0 og x1 eru tolur sem f raedur vid
#           n er tala
# Gildi:    xn er rot n-ta talan
#           i runu sem nalgar rot f,
#           reiknud með snidilsadferd
def secant_method(f,x0,x1,n):
    if n == 1:
        return x1
    else:
        x2 = x1 - (f(x1)*(x1-x0))/(f(x1)-f(x0))
        return secant_method(f,x1,x2,n-1)


# Notkun:   false_position(f,a,b)
# Fyrir:    f er samfelld einundaradgerd
#           skilgreind á  bilinu [a,b].
#           a < b eru tölur þannig að f(a)*f(b) < 0
#           n er jakvaed heiltala
# Eftir:    (c,fc) tuple þar sem c er rot f og fc er f(c)
def false_position(f,a,b,n):
    fa=f(a)
    fb=f(b)
    for i in range (n,0,-1):
        c=(b*fa - a*fb)/(fa-fb)
        fc=f(c)
#        print(fc)
        if fa*fc <= 0:
            b=c
            fb=fc
        else:
            a=c
            fa=fc
    return c,fc
    
# andhverf kvadratisk bruun
def IQI(f,x0,x1,x2,n):
    if n == 0:
        return x2
    else:
        fx0, fx1, fx2 = f(x0), f(x1), f(x2)
        q = fx0/fx1
        r = fx2/fx1
        s = fx2/fx0
        num = r*(r-q)*(x2-x1) + (1-r)*s*(x2-x0)
        den = (q-1)*(r-1)*(s-1)
        x3 = x2 - num/den
        return IQI(f,x1,x2,x3,n-1)
    
        
x0, x1, x2 = 1, 2, 0

f = lambda x: math.exp(x) + math.sin(x) - 4
x2sm = secant_method(f,x0,x1,2)
x3sm = secant_method(f,x0,x1,3)

x2fp = false_position(f,x0,x1,2)[0]

x2iqi = IQI(f,x0,x1,x2,2)

print("rot f(x) m.v. snidilsadferd og eina itrun er:                  ", x2sm)
print("rot f(x) m.v. snidilsadferd og tvaer itranir er:               ", x3sm)
print("rot f(x) m.v. adferd rangrar stodu og tvaer itranir er:        ", x2fp)
print("rot f(x) m.v. andhverfa kvadratiska bruun og tvaer itranir er: ", x2iqi)
# rot f(x) m.v. snidilsadferd og eina itrun er:                   1.0929065801160904
# rot f(x) m.v. snidilsadferd og tvaer itranir er:                1.1193566855644101
# rot f(x) m.v. adferd rangrar stodu og tvaer itranir er:         1.1193566855644101
# rot f(x) m.v. andhverfa kvadratiska bruun og tvaer itranir er:  1.1292724601823607


##### Dæmi 2
def fpitr(g,x,n):
    for i in range(n):
        x = g(x)
        print(x)
        
r = 2**(1/2)
g1 = lambda x: (x/2 + 1/x)
g2 = lambda x: (2*x/3 + 2/(3*x))
g3 = lambda x: (3*x/4 + 1/(2*x))
g4 = lambda x: np.cos(x)

lausn221 = fpitr(g1,r,10)


# Dæmi 1
# Notkun:   helmadf(f,a,b)
# Fyrir:    f er samfellt fall skilgreint á  bilinu [a,b]
#           a < b eru tölur þannig að f(a)*f(b) < 0
#           TOL er tala
# Eftir:    (c,fc) tuple þar sem c er rót f og fc er f(c)
def helmadf(f,a,b,TOL):
    fa=f(a)
    fb=f(b)
    while (b-a)/2 > TOL:
        c=(a+b)/2
        fc=f(c)
#        print(fc)
        if fa*fc <= 0:
            b=c
            fb=fc
        else:
            a=c
            fa=fc
    return c,fc
    

# Skilgreinum fallið sem við viljum finna rót á
f = lambda x: (np.linalg.det((np.array([[1,2,3,x],[4,5,x,6],[7,x,8,9],[x,10,11,12]]))) - 1000)
# Prófum ýmsar prentanir til að finna möguleg gildi a og b í helmingunaraðferðinni, m.a.
# for i in range(-20,20,2):
#     print(i, " ... ", f(i))
# Veljum: 
a1 = -18
b1 = -16
a2 = 8
b2 = 10
TOL = 0.5e-6
lausn1 = helmadf(f,a1,b1,TOL)
# Out[30]: (-17.188498497009277, 0.004138882266488508)
print("x1 = %.6f" % lausn1[0])
# x1 = -17.188498
lausn2 = helmadf(f,a2,b2,TOL)
# Out[32]: (9.708298683166504, -0.0005025448032256463)
print("x2 = %.6f" % lausn2[0])
# x2 = 9.708299

# Exercise 1.1 Computer Problems, Dæmi 9
ex11c9 = lambda x: (np.pi)*x*x*(1-x/3) - 1
helmadf(ex11c9,0,1,1e-4)


########## Dæmablað 1 ##########
# Dæmi 1
# Notkun:   nested_multiplication(d,c,x,b)
# Fyrir:    d er heiltala, stig margliðunnar
#           c er tuple, stuðlar margliðunnar
#           x er tala
#           b er valfjáls tuple
# Eftir:    gildi margliðunnar í punktinum x
def nested_multiplication(d,c,x,b=0):
    if b == 0:
        b=np.zeros(d+1)
    y1 = c[d]
    for i in range(d,0,-1):
        y1 = c[i]+y1*(x-b[i])
    return y1
    
# Búum til tuple með 51 ás
a = ()

for i in range(0,51):
    a += (1,)
    
# Reiknum 1.00001^0 + 1.00001^1 + ... + 1.00001^50
P = nested_multiplication(50,a,1.00001)


# Búum til tuple með 51 ás
b = ()

for i in range(0,100):
    if i % 2 == 0:
        b += (1,)
    else:
        b += (-1,)

# Reiknum 1.00001^0 - 1.00001^1 + 1.00001^2 - ... + 1.00001^98 - 1.00001^99
Q = nested_multiplication(97,b,1.00001)

# Prentum niðurstöður
print(P)
# 51.01275208274999

print(Q)
# -2.001480725434141


### Dæmi 2
a = -12345678987654321
b = 123

ar_minus = -111111111
ar_plus  =  111111111

c = (ar_plus * ar_minus)
print(a==c)
#True

#Sé að a og c eru geymdar sem -12345678987654320 svo ég prófa þetta... áhugavert
d = -12345678987654320
print(a==d)
#False

#Lausn verkefnis:
ans1 = a + (a*a + b*b)**(1/2)

print("Verri aðferð prentað sem f með 2 aukastöfum: %.4f" % (ans1))
# Verri aðferð prentað sem f með 2 aukastöfum: 0.0000

print("Verri aðferð prentað sem e með 2 aukastöfum: %.4e" % (ans1))
# Verri aðferð prentað sem e með 2 aukastöfum: 0.0000e+00

ans2 = -b*b / (a-(a*a + b*b)**(1/2))

print("Betri aðferð prentað sem f með 2 aukastöfum: %.4f" % (ans2))
# Betri aðferð prentað sem f með 2 aukastöfum: 0.0000

print("Betri aðferð prentað sem e með 2 aukastöfum: %.4e" % (ans2))
# Betri aðferð prentað sem e með 2 aukastöfum: 6.1272e-13
