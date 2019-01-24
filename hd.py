import numpy as np 
import math

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
