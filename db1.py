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