import numpy as np 

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
