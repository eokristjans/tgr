a = 2 
b = 3

def fun(x, a = a, b = b):
    return x + a + b

# Hér ætti að koma út 1 + 2 + 3 = 6
print(fun(1))

a = 0 
b = 0

# Hér ætti að koma út 1 + 0 + 0 = 1
# en það kemur út 6
print(fun(1))

################################

print("\n Gerum þetta með því að nota global lykilorðið \n")

global c, d
c = 2
d = 3 

def fun2(y, c = c, d = d):
    return y + a + b

# Hér ætti að koma út 1 + 2 + 3 = 6
print(fun(1))

c = 0 
d = 0

# Hér ætti að koma út 1 + 0 + 0 = 1
print(fun2(1))

print("\n Sleppum að nota default arguments \n")

y = 1
e = 2
f = 3 

def fun3(y, e, f):
    return y + e + f

# Hér ætti að koma út 1 + 2 + 3 = 6
print(fun3(y, e, f))

e = 0 
f = 0

# Hér ætti að koma út 1 + 0 + 0 = 1
print(fun3(y, e, f))